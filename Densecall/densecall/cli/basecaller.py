import torch
import numpy as np
import os, sys, re, argparse, pickle, time, glob
from torch.utils.data import Dataset, DataLoader
from fast_ctc_decode import beam_search
from ont_fast5_api.fast5_interface import get_fast5_file
from torch.multiprocessing import Queue, Process
from densecall.ctc import Model
from collections import OrderedDict
import toml
from argparse import ArgumentDefaultsHelpFormatter
from densecall.util import __models__

def med_mad(x, factor=1.4826):
    med = np.median(x)
    mad = np.median(np.absolute(x - med)) * factor
    return med, mad


def convert_statedict(state_dict):
    from collections import OrderedDict
    new_checkpoint = OrderedDict()
    for k, v in state_dict.items():
        name = k[7:] if k.startswith('module.') else k # remove module.
        new_checkpoint[name] = v
    return new_checkpoint


def load_model(args = None):
    dirname = args.model_directory
    if not os.path.isdir(dirname) and os.path.isdir(os.path.join(__models__, dirname)):
        dirname = os.path.join(__models__, dirname)
    modelfile = find_mat_recursive(dirname, ends='.tar')[0]
    config = toml.load(os.path.join(dirname, 'config.toml'))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    network = Model(config=config)
    state_dict = torch.load(modelfile, map_location=device)
    network.load_state_dict(convert_statedict(state_dict))
    if args.amp:
        network.half() 

    network = network.to(device)
    network.eval()
    torch.set_grad_enabled(False)
    return network, device


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


def trim(signal, window_size=40, threshold_factor=2.4, min_elements=3):
    """

    from: https://github.com/nanoporetech/bonito/blob/master/bonito/fast5.py
    """

    min_trim = 10
    signal = signal[min_trim:]

    med, mad = med_mad(signal[-(window_size*100):])

    threshold = med + mad * threshold_factor
    num_windows = len(signal) // window_size

    seen_peak = False

    for pos in range(num_windows):
        start = pos * window_size
        end = start + window_size
        window = signal[start:end]
        if len(window[window > threshold]) > min_elements or seen_peak:
            seen_peak = True
            if window[-1] > threshold:
                continue
            return min(end + min_trim, len(signal)), len(signal)

    return min_trim, len(signal)

def load_read_data(read_file, platform='kq'):
    """
    #awesomely here,use normalize=True,I am not sure why it does not use median normalization

    Args:
        read_file ([type]): [description]
        platform (str, optional): [description]. Defaults to 'kq'.

    Returns:
        [type]: [description]
    """
    f5 = get_fast5_file(read_file, mode="r")
    for read in f5.get_reads():
        try:
            signal = read.get_raw_data(scale=True)

        except:
            continue
        signal_start = 0
        signal_end = len(signal)
        med, mad = med_mad(signal[signal_start:signal_end])
        norm_signal = (signal[signal_start:signal_end] - med) / mad
        yield norm_signal, read.read_id




def basecall(read_file_name, fo, args, network, device):
    basename = os.path.basename(read_file_name)
    try:
        for norm_signal, read_id in load_read_data(read_file_name):
            fa, qstring = predict(norm_signal, args, network, device)
            print(f'@{read_id}\n{fa}\n+\n{qstring}', file=fo)
    except Exception as e:
        print(e)
        print("error at file", read_file_name)
        return


def segment_with_overlap(seg, s, overlap=0):
    assert overlap < s
    seg = np.concatenate((seg, np.zeros((s - overlap))))
    nrows = (seg.size - s) // (s - overlap) + 1
    n = seg.strides[0]
    return np.lib.stride_tricks.as_strided(seg, shape=(nrows, s), strides=((s - overlap) * n, n))


def concat(xs, dim=0):
    """
    Type agnostic concat.
    """
    if isinstance(xs[0], torch.Tensor):
        return torch.cat(xs, dim=dim)
    elif isinstance(xs[0], np.ndarray):
        return np.concatenate(xs, axis=dim)
    elif isinstance(xs[0], list):
        return [x for l in xs for x in l]
    elif isinstance(xs[0], str):
        return ''.join(xs)
    elif isinstance(xs[0], dict):
        return {k: concat([x[k] for x in xs], dim) for k in xs[0].keys()}
    else:
        raise TypeError
    
    
def stitch(chunks, chunksize, overlap, length, stride, reverse=False):
    """
    Stitch chunks together with a given overlap
    Copied from: https://github.com/nanoporetech/bonito/blob/master/bonito/cli/basecaller.py
    """
    if chunks.shape[0] == 1: return chunks.squeeze(0)

    semi_overlap = overlap // 2
    start, end = semi_overlap // stride, (chunksize - semi_overlap) // stride
    stub = (length - overlap) % (chunksize - overlap)
    first_chunk_end = (stub + semi_overlap) // stride if (stub > 0) else end

    if reverse:
        chunks = list(chunks)
        return concat([
            chunks[-1][:-start], *(x[-end:-start] for x in reversed(chunks[1:-1])), chunks[0][-first_chunk_end:]
        ])
    else:
        return concat([
            chunks[0, :first_chunk_end], *chunks[1:-1, start:end], chunks[-1, start:]
        ])


def permute(x, input_layout, output_layout):
    """
    Permute `x` from `input_layout` to `output_layout`

    >>> permute(x, 'TNC', 'NTC')
    """
    if input_layout == output_layout: return x
    return x.permute(*[input_layout.index(x) for x in output_layout])


def predict(signal, args, network, device):

    chunks = segment_with_overlap(signal, args.seqlen, args.overlap)
    scores = []
    for i in range(0, chunks.shape[0], args.batchsize):
        end = i+args.batchsize
        if end > chunks.shape[0]: end = chunks.shape[0]
        event = torch.unsqueeze(torch.FloatTensor(chunks[i:end]), 1)
        event = event.to(device, non_blocking=True) 
        out=network.forward(event.half() if args.amp else event) # (T, N, C)
        out = permute(out, 'TNC', 'NTC')
        scores.append(out.cpu().to(torch.float32))

    scores = concat(scores)
    scores = stitch(scores, args.seqlen, args.overlap, len(signal), network.stride) 
    seq = network.decode(scores, beamsize=5)
    qstring = '*'*len(seq)
    if args.reverse:
        seq = seq[::-1]
        qstring = qstring[::-1]
    return seq, qstring

    
def real_main(args, network, device):
    fo = open(args.output_file, 'w')
    all_fast5 = find_mat_recursive(args.fast5dir)
    
    poem = np.arange(len(all_fast5))
    np.random.seed(1)
    np.random.shuffle(poem)
    all_fast5 = [all_fast5[i] for i in poem][:args.max_reads]
    
    for fast5 in all_fast5:
        basecall(fast5, fo, args, network, device)
    fo.close()

def argparser():
    parser = argparse.ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False)
    
    parser.add_argument("model_directory")
    parser.add_argument("fast5dir", default=None, type=str)
    parser.add_argument("--output-file", type=str, help='output fastq file')
    parser.add_argument("-r", "--reverse", default=True, action="store_true", help="reverse for RNA (default: True)")
    parser.add_argument("-seqlen", "--seqlen", default=4096, type=int, help="chunksize (default: 4096)")
    parser.add_argument("-overlap", "--overlap", default=200, type=int, help="overlap between chunks (default: 200)")
    parser.add_argument("-b", "--batchsize", default=16, type=int, help="default: 16")
    parser.add_argument("-B", "--beam-size", default=5, type=int, help="CTC beam search size (default: 5)")
    parser.add_argument("-amp", "--amp", default=False, action="store_true")
    parser.add_argument("--max-reads", type=int, default=None, help='Number of reads to use (default: all)')
    
    return parser


def main(args):
    torch.multiprocessing.set_start_method('spawn')
    
    torch.backends.cudnn.enabled = True
    torch.backends.cudnn.deterministic = True
    network, device = load_model(args)
    real_main(args, network, device)
    

    
