import h5py
import os
import re
import torch
from torch.utils.data import Dataset, DataLoader
import parasail, re
from collections import defaultdict
from torch.cuda.amp import autocast, GradScaler
from pytorch_ranger import Ranger
import numpy as np
from densecall.ctc import model 
import tempfile
import warnings
import time
import importlib
import argparse
from argparse import ArgumentDefaultsHelpFormatter
import toml


warnings.filterwarnings("ignore", category=UserWarning)



    
split_cigar = re.compile(r"(?P<len>\d+)(?P<op>\D+)")


def find_mat_recursive(wk, ends='.fast5'):
    """
    wk:directory, which contain mat files\n
    auto find all file ends with .mat
    """
    mat_lst = []
    for root, filefolder, filelist in os.walk(wk):
        for mat in filelist:
            if mat.lower().endswith(ends):
                filename = os.path.join(root, mat)
                mat_lst.append(filename)
    return mat_lst

def read_signal(signal):
    value_m = np.median(signal)
    mad = np.median(np.abs(signal - value_m)) * 1.4826
    signal = (signal - value_m) / mad
    return signal

def parasail_to_sam(result, seq):
    """
    Extract reference start and sam compatible cigar string.

    :param result: parasail alignment result.
    :param seq: query sequence.

    :returns: reference start coordinate, cigar string.
    """
    cigstr = result.cigar.decode.decode()
    first = re.search(split_cigar, cigstr)

    first_count, first_op = first.groups()
    prefix = first.group()
    rstart = result.cigar.beg_ref
    cliplen = result.cigar.beg_query

    clip = '' if cliplen == 0 else '{}S'.format(cliplen)
    if first_op == 'I':
        pre = '{}S'.format(int(first_count) + cliplen)
    elif first_op == 'D':
        pre = clip
        rstart = int(first_count)
    else:
        pre = '{}{}'.format(clip, prefix)

    mid = cigstr[len(prefix):]
    end_clip = len(seq) - result.end_query - 1
    suf = '{}S'.format(end_clip) if end_clip > 0 else ''
    new_cigstr = ''.join((pre, mid, suf))
    return rstart, new_cigstr


def accuracy(ref, seq, balanced=False, min_coverage=0.0):
    """
    Calculate the accuracy between `ref` and `seq`
    """
    alignment = parasail.sw_trace_striped_32(seq, ref, 8, 4, parasail.dnafull)
    counts = defaultdict(int)

    q_coverage = len(alignment.traceback.query) / len(seq)
    r_coverage = len(alignment.traceback.ref) / len(ref)

    if r_coverage < min_coverage:
        return 0.0

    _, cigar = parasail_to_sam(alignment, seq)

    for count, op in re.findall(split_cigar, cigar):
        counts[op] += int(count)

    if balanced:
        accuracy = (counts['='] - counts['I']) / (counts['='] + counts['X'] +
                                                  counts['D'])
    else:
        accuracy = counts['='] / (counts['='] + counts['I'] + counts['X'] +
                                  counts['D'])
    return accuracy * 100


def index2base(read, alphabet='NACGT'):
    """Transfer the number into dna base.
    The transfer will go through each element of the input int vector.
    """

    bpread = [alphabet[x] for x in read]
    bpread = ''.join(x for x in bpread)
    return bpread



def validate(bpreads, y_true, y_len):
    reference = [
                index2base(y[:y_len[inx]])
                for inx, y in enumerate(y_true)
            ]

    accuracy_lst = [
        accuracy(ref, seq, min_coverage=0.5) if len(seq) else 0.
        for ref, seq in zip(reference, bpreads)
    ]
    return accuracy_lst

def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)
    
class dataloader(Dataset):
    def __init__(self, *data, index=False):
        self.event, self.event_len, self.label, self.label_len = data
        self.len = len(self.event)
        self.index = index

    def __getitem__(self, index):
        event = self.event[index]
        event_len = self.event_len
        label = self.label[index]
        label_len = self.label_len[index]
        if not self.index:
            return (event, event_len, label, label_len)
        else:
            return (event, event_len, label, label_len, index)

    def __len__(self):
        return self.len   

def save(**metrics):
        """Save the model state
        """

        torch.save(metrics['network'].state_dict(), os.path.join(metrics['workdir'], f"weights_{metrics['epoch']}_{metrics['accu']:.2f}.tar"))
        #torch.save({'optimizer_state': metrics['optimizer'].state_dict()}, os.path.join(metrics['workdir'], f"optim_{metrics['epoch']}.tar"))
        #torch.save({'lr_scheduler_state':metrics['scheduler']}, os.path.join(metrics['workdir'], f"scheduler_{metrics['epoch']}.tar"))
        toml.dump(metrics['config'], open(os.path.join(metrics['workdir'], 'config.toml'), 'w'))
        


def load_weights(workdir, checkpoint_file=None):
    if checkpoint_file:
        weights_file = os.path.join(workdir, checkpoint_file)
        start_epoch = os.path.basename(weights_file).split('_')[1]
        print('loading weights:', weights_file)
        return weights_file, int(start_epoch)
    weights_files = find_mat_recursive(workdir, ends='.tar')
    weights_files = [f for f in weights_files if re.search(r'weights.*\.tar$', f)]
    weights_file = sorted(weights_files, key=lambda s: int(os.path.basename(s).split('_')[1]))[-1]
    start_epoch = os.path.basename(weights_file).split('_')[1]
    
    print('loading weights:', weights_file)
    return weights_file, int(start_epoch)
    
def load_data(data_dir):
    
    x_train = np.load(f'{data_dir}/chunks.npy')
    y_train = np.load(f'{data_dir}/references.npy').astype('int64')
    x_len_train = np.full((x_train.shape[0], ), x_train.shape[1], dtype='int64')
    y_len_train = np.sum(y_train != 0, axis = 1).astype('int64')
    
    np.random.seed(42)
    indices = np.random.permutation(len(x_train))
    x_train = x_train[indices]
    y_train = y_train[indices]
    x_len_train = x_len_train[indices]
    y_len_train = y_len_train[indices]

    print("train dataset:", x_train.shape, y_train.shape, x_len_train.shape, y_len_train.shape)
    
    if os.path.exists(os.path.join(data_dir, 'validation')):
        print('loading validation')
        
        x_valid = np.load(f'{data_dir}/validation/chunks.npy')
        y_valid = np.load(f'{data_dir}/validation/references.npy').astype('int64')
        x_len_valid = np.full((x_valid.shape[0], ), x_valid.shape[1], dtype='int64')
        y_len_valid = np.sum(y_valid != 0, axis = 1).astype('int64')
    else:
        print('splitting from train dataset into validation')
        
        n = int(len(x_train)*0.9)
        x_train, x_valid = x_train[:n], x_train[n:]
        y_train, y_valid = y_train[:n], y_train[n:]
        x_len_train, x_len_valid = x_len_train[:n], x_len_train[n:]
        y_len_train, y_len_valid = y_len_train[:n], y_len_train[n:]

    print("valid dataset:", x_valid.shape, y_valid.shape, x_len_valid.shape, y_len_valid.shape)
    
    return x_train, x_valid, y_train, y_valid, x_len_train, x_len_valid, y_len_train, y_len_valid


class dataloader_hdf5(Dataset):
    def __init__(self, recfile="/tmp/train.hdf5", seq_len=4096, index=False, elen=342):
        self.recfile = recfile
        self.seq_len = seq_len
        self.index = index
        h5 = h5py.File(self.recfile, "r")
        self.len = len(h5["events"])
        h5.close()
        self.elen = elen
        print("Dataloader total events:", self.len, "seqlen:", self.seq_len, "event len:", self.elen)

    def __getitem__(self, index):
        h5 = h5py.File(self.recfile, "r")
        event = h5["events"][index]
        event_len = self.elen
        label = h5["labels"][index]
        label_len = h5["labels_len"][index]
        h5.close()
        if not self.index:
            return (event, event_len, label, label_len)
        else:
            return (event, event_len, label, label_len, index)

    def __len__(self):
        return self.len
    
def convert_statedict(state_dict):
    from collections import OrderedDict
    new_state_dict = OrderedDict()
    for k, v in state_dict.items():
        if k.startswith('module'):
            name = k[7:]
            new_state_dict[name] = v
        else:
            new_state_dict[k] = v
    return new_state_dict  

def train(config, args):
    
    os.makedirs(args.workdir, exist_ok=True)
    
    output_classes = 5  

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    start_epoch = 0
    scaler = GradScaler()

    if not args.retrain:
        network = model.Model(config=config).to(device)
    else:
        print('retrain network, load prev config')
        config = toml.load(os.path.join(args.workdir, 'config.toml'))
        
        network = model.Model(config=config)
        network = network.to(device) 
        start_epoch = args.start_epoch
        state_dict = torch.load(args.checkpointfile, map_location=device)
        network.load_state_dict(convert_statedict(state_dict))
        
    print(config)
    if args.fine_tuning: 
        print('fine tuning,frozen all layers paramerter other than decoder')
        for name, param in network.named_parameters():
            if 'decoder' not in name: # decoder 将被微调
                param.requires_grad = False
        network = network.to(device) 

    #network = torch.compile(network)
    
    
    network.eval()
    with torch.no_grad():
        fakedata = torch.randn(args.batch_size, 1, args.seqlen) 
        fakeout = network.forward(fakedata.to(device))
        event_len = fakeout.shape[0]
        print('output shape:', fakeout.shape) 
     
    if args.npy_data:   
        x_train, x_valid, y_train, y_valid, x_len_train, x_len_valid, y_len_train, y_len_valid = load_data(args.data_dir)
        data_train = dataloader(x_train, event_len, y_train, y_len_train)
        data_valid = dataloader(x_valid, event_len, y_valid, y_len_valid)
        
    else:
        data_train = dataloader_hdf5(recfile=os.path.join(args.data_dir,'train.hdf5'), seq_len=args.seqlen, elen=event_len)
        data_valid = dataloader_hdf5(recfile=os.path.join(args.data_dir,'valid.hdf5'), seq_len=args.seqlen, elen=event_len)
        
    data_loader_train = DataLoader(dataset=data_train, batch_size=args.batch_size, shuffle=True, num_workers=4, pin_memory=True)
    data_loader_valid = DataLoader(dataset=data_valid, batch_size=args.batch_size, shuffle=True, num_workers=4, pin_memory=True)  

    num_params = count_parameters(network)
    print(f'The network has {num_params} trainable parameters.')

    optimizer = Ranger(filter(lambda p: p.requires_grad, network.parameters()), lr=args.lr, weight_decay=0.01) 
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=1, factor=0.5, threshold=0.1, min_lr=1e-05)   
    
    smoothweights = torch.cat([torch.tensor([0.1]), (0.1 / (output_classes - 1)) * torch.ones(output_classes - 1)]).to(device)  

    
    
    count = 0  
    torch.cuda.empty_cache()
    print('total parameters:', sum(p.numel() for p in network.parameters() if p.requires_grad))
    for epoch in range(start_epoch, args.epochs):
        network.train()
        loopcount = 0
        totalloss = 0
        start2 = time.time()
        for i, (event, event_len, targets, target_lengths) in enumerate(data_loader_train): 
            start = time.time()
            event = torch.unsqueeze(event, 1)
            if len(event) != args.batch_size: continue
            targets = targets[:, :max(target_lengths)]
            event = event.float().to(device, non_blocking=True)
            targets = targets.to(device, non_blocking=True)
            event_len = event_len.to(device, non_blocking=True)
            target_lengths = target_lengths.to(device, non_blocking=True)

            optimizer.zero_grad()

            with autocast():
                log_probs = network.forward(event)
                loss, _ = network.calculate_loss(targets, log_probs.to(torch.float32))

            label_smoothing_loss = -((log_probs * smoothweights.to(log_probs.device)).mean())
            label_smoothing_loss += loss
            scaler.scale(label_smoothing_loss).backward()

            end = time.time()
            speed = int(args.batch_size / (end - start))
            totalloss+=loss.cpu().detach().numpy()
            if i % 100 == 0: print("Loss:", loss.data, "epoch:", epoch, count, "chunks/s:", speed, "lr:", optimizer.param_groups[0]['lr'])
            
            scaler.unscale_(optimizer)
            grad_norm = torch.nn.utils.clip_grad_norm_(network.parameters(), 2.0, error_if_nonfinite=False)
            scaler.step(optimizer)
            scaler.update()
            loopcount +=1
            count +=1

            if loopcount >= args.train_loopcount: break
            #torch.cuda.empty_cache()
        end2 = time.time()
        duration = end2 - start2
        train_mean_loss = totalloss/loopcount
        print("Train epoch loss", train_mean_loss, f"Train epoch cost {duration:2f} s")
        

        val_total  = 0
        valid_totalloss = 0
        accuracy_lst = []
        seqs_lst = []
        network.eval() 
        with torch.no_grad():
            
            for i, (val_event, val_event_len, val_targets, val_target_lengths) in enumerate(data_loader_valid): 
                val_event = torch.unsqueeze(val_event, 1)
                if len(val_event) != args.batch_size: continue
                val_event = val_event.float().to(device, non_blocking=True)
                val_targets = val_targets.to(device, non_blocking=True)
                with autocast():
                    val_log_probs = network.forward(val_event)
                    val_loss, _ = network.calculate_loss(val_targets, val_log_probs.to(torch.float32))
                    seqs = network.basecall(val_log_probs.to(torch.float32), greedy=True)
                    seqs_lst += seqs
                    accuracy_lst += validate(seqs, val_targets, val_target_lengths)
                    
                valid_totalloss += val_loss.cpu().detach().numpy()
                val_total +=1
                if val_total > args.batch_size:break
        valid_accu = np.mean(accuracy_lst)
        valid_mean_loss = valid_totalloss / val_total
        scheduler.step(valid_mean_loss)    
        print("valid epoch accuracy", valid_accu, "loss", valid_mean_loss)
        

        config['optimizer']['lr'] = args.lr
        metrics = {'loss':train_mean_loss, 'valid_loss':valid_mean_loss, 'accu':valid_accu,
                   'network':network, 'scaler':scaler, 'optimizer':optimizer, 'scheduler':scheduler,
                   'workdir':args.workdir, 'config':config,'epoch':epoch,
                   
                   }
        save(**metrics)
 
        
        
def argparser():
    parser = argparse.ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False)

    parser.add_argument('--config', type=str, default='rna-config.toml', 
                        help='Path to the configuration file (default: rna-config.toml)')
    parser.add_argument('--workdir', type=str, default='model', 
                        help='Directory to save the trained model (default: model)')
    parser.add_argument('--data_dir', type=str, default='../data/rna', 
                        help='Directory containing the training data, which contains train.hdf5 and valid.hdf5')
    parser.add_argument('--checkpointfile', type=str, default=None, 
                        help='Path to the checkpoint file to load model from (default: None)')
    parser.add_argument('--start_epoch', type=int, default=0, 
                        help='start new epoch from this when fine tuning model (default: 0)')
    parser.add_argument('--retrain', action='store_true', 
                        help='Whether to retrain the model from scratch (default: False)')
    parser.add_argument('--batch_size', type=int, default=32, 
                        help='Batch size for training (default: 32)')
    parser.add_argument('--lr', type=float, default=0.002, 
                        help='Learning rate for training (default: 0.002)')
    parser.add_argument('--epochs', type=int, default=30, 
                        help='Number of epochs for training (default: 30)')
    parser.add_argument('--train_loopcount', type=int, default=1000000, 
                        help='Number of training loop counts (default: 1000000)')
    parser.add_argument('--fine_tuning', default=False, action='store_true',
                        help='fine tune one of the pretrained models (default: False)')
    parser.add_argument('--npy_data', default=False, action='store_true')
    

    return parser

def main(args):
    config = toml.load(args.config)
    torch.backends.cudnn.benchmark = True
    args.seqlen = config['basecaller']['chunksize']
    train(config, args)        
        

