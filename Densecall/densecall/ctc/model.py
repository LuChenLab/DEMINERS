import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np 
from fast_ctc_decode import beam_search, viterbi_search, crf_greedy_search, crf_beam_search
import sys
import torch.utils.checkpoint as cp
from collections import OrderedDict


class Mish(nn.Module):
    def __init__(self):
        super().__init__()
    def forward(self, x):
        return x *( torch.tanh(torch.nn.functional.softplus(x)))

class squeeze_excite(torch.nn.Module):
    def __init__(self, in_channels = 512, size=1, reduction="/16", activation=torch.nn.GELU):
        super(squeeze_excite, self).__init__()
        self.in_channels = in_channels
        self.avg = torch.nn.AdaptiveAvgPool1d(1)
        if type(reduction) == str:
            self.reductionsize = self.in_channels // int(reduction[1:])
        else:
            self.reductionsize = reduction
        self.fc1 = nn.Linear(self.in_channels, self.reductionsize)
        self.activation = activation() # was nn.ReLU(inplace=True)
        self.fc2 = nn.Linear(self.reductionsize, self.in_channels)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        input = x
        x = self.avg(x)
        x = x.permute(0,2,1)
        x = self.activation(self.fc1(x))
        x = self.sigmoid(self.fc2(x))
        return input * x.permute(0,2,1)



    
class TemporalBlock(nn.Module):
    def __init__(self, n_inputs, n_outputs, kernel_size, stride, dilation, padding, dropout=0.2):
        super(TemporalBlock, self).__init__()
        self.conv1 = nn.Conv1d(in_channels=n_inputs,
                               out_channels=n_outputs,
                               kernel_size=kernel_size,
                               stride=stride,
                               padding=padding,
                               dilation=dilation)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(dropout)
        
        self.net = nn.Sequential(
            self.conv1,
            self.relu,
            self.dropout
        )
        self.downsample = nn.Conv1d(n_inputs, n_outputs, 1) if n_inputs != n_outputs else None

    def forward(self, x):
        out = self.net(x)

        res = x if self.downsample is None else self.downsample(x)
        return self.relu(out + res)

class TemporalConvNet(nn.Module):
    def __init__(self, num_inputs, num_channels, kernel_size=3, dropout=0.2):
        super(TemporalConvNet, self).__init__()
        layers = []
        num_levels = len(num_channels)
        for i in range(num_levels):
            dilation_size = 2 ** i
            print(dilation_size)
            in_channels = num_inputs if i == 0 else num_channels[i-1]
            out_channels = num_channels[i]
            layers += [TemporalBlock(in_channels, out_channels, kernel_size, stride=1, dilation=dilation_size,
                                      padding=(kernel_size-1) * dilation_size // 2, dropout=dropout)]

        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return self.network(x)

class BiTemporalConvNet(nn.Module):
    def __init__(self, num_inputs, num_channels, kernel_size=2, dropout=0.2):
        super(BiTemporalConvNet, self).__init__()
        self.forward_tcn = TemporalConvNet(num_inputs, num_channels, kernel_size, dropout)
        self.backward_tcn = TemporalConvNet(num_inputs, num_channels, kernel_size, dropout)

    def forward(self, x):
        x_reversed = torch.flip(x, dims=[2])
        out_forward = self.forward_tcn(x)
   
        out_backward = self.backward_tcn(x_reversed)
        out_backward = torch.flip(out_backward, dims=[2])
        
        # Concatenate the forward and backward outputs
        out = torch.cat((out_forward, out_backward), dim=1)
        return out
    
    
def _bn_function_factory(norm, relu, conv):
    def bn_function(*inputs):
        concated_features = torch.cat(inputs, 1)
        bottleneck_output = conv(relu(norm(concated_features)))
        return bottleneck_output

    return bn_function


class _DenseLayer(nn.Module):
    def __init__(self, num_input_features, growth_rate, bn_size, drop_rate, kernel_size=3, activation=nn.SiLU):
        super(_DenseLayer, self).__init__()

        self.add_module('norm1', nn.BatchNorm1d(num_input_features))
        self.add_module('act1', activation())
        self.add_module('conv1', nn.Conv1d(num_input_features, bn_size * growth_rate,
                                        kernel_size=1, stride=1, bias=False))
        self.add_module('norm2', nn.BatchNorm1d( bn_size * growth_rate))
        self.add_module('act2', activation())
        self.add_module('conv2', nn.Conv1d(bn_size * growth_rate, growth_rate,
                                        kernel_size=kernel_size, stride=1, padding=kernel_size//2, bias=False))

        self.drop_rate = drop_rate
        self.efficient = True

    def forward(self, *prev_features):
        
        bn_function = _bn_function_factory(self.norm1, self.act1, self.conv1)
        if self.efficient and any(prev_feature.requires_grad for prev_feature in prev_features):
            bottleneck_output = cp.checkpoint(bn_function, *prev_features)
        else:
            bottleneck_output = bn_function(*prev_features)
        new_features = self.conv2(self.act2(self.norm2(bottleneck_output)))
        if self.drop_rate > 0:
            new_features = F.dropout(new_features, p=self.drop_rate, training=self.training)
        return new_features


    
class _DenseBlock(nn.Module):
    def __init__(self, num_layers, num_input_features, bn_size, growth_rate, drop_rate,kernel_size=3, kernel_size_increment=0, activation=nn.SiLU):
        super(_DenseBlock, self).__init__()

        for i in range(num_layers):
            layer_kernel_size = kernel_size + i * kernel_size_increment 
            if layer_kernel_size % 2 == 0:
                layer_kernel_size += 1
            if layer_kernel_size > 100:
                layer_kernel_size = 99

            layer=_DenseLayer(num_input_features + i * growth_rate, growth_rate, bn_size, drop_rate,kernel_size=layer_kernel_size, activation=activation)
            self.add_module('denselayer%d' % (i + 1), layer)
        
    def forward(self, init_features):
        features = [init_features]
        for name, layer in self.named_children():
            new_features = layer(*features)
            features.append(new_features)
        return torch.cat(features, 1)


class _Transition(nn.Sequential):
    def __init__(self, num_input_features, num_output_features, stride=2, activation=nn.SiLU):
        super(_Transition, self).__init__()
        self.add_module('norm', nn.BatchNorm1d(num_input_features))
        self.add_module('act', activation())
        self.add_module('conv', nn.Conv1d(num_input_features, num_output_features,
                                          kernel_size=1, stride=1, bias=False))
        if stride != 1:
            self.add_module('pool', nn.AvgPool1d(kernel_size=stride, stride=stride, padding=stride//2))

def gcd(a, b):
    while b != 0:
        temp = b
        b = a % b
        a = temp
    return a
    
    
class StackedBidirectionalGRU(nn.Module):
    def __init__(self, input_size, hidden_sizes, num_layers):
        super(StackedBidirectionalGRU, self).__init__()
        

        if len(hidden_sizes) != num_layers:
            raise ValueError("Length of hidden_sizes must be equal to num_layers")


        self.stacked_gru = nn.ModuleList()
        for i in range(num_layers):
            input_dim = input_size if i == 0 else 2 * hidden_sizes[i-1]  
            self.stacked_gru.append(nn.GRU(input_size=input_dim,
                                           hidden_size=hidden_sizes[i],
                                           num_layers=1,  
                                           bidirectional=True,
                                           batch_first=True))  

    def forward(self, x):
        for i, gru_layer in enumerate(self.stacked_gru):
            x, _ = gru_layer(x)
        return x
    
    
class Model(nn.Module):
    def __init__(self, config={}, activation=nn.SiLU):

        # load hyper param
        self.config = config
        
        # decoder
        
        # important
        super(Model, self).__init__()
        
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.criterions = {}
        self.criterions['ctc'] = nn.CTCLoss(blank = 0, zero_infinity = True).to(self.device)
        # cnn
        self.growth_rate = config['cnn']['growth_rate']
        self.block_config = config['cnn']['block_config']
        self.num_init_features = config['cnn']['num_init_features']
        self.bn_size = config['cnn']['bn_size']
        self.drop_rate = config['drop']['drop_rate']

        self.kernel_size = config['cnn']['kernel_size']
        self.stride = config['cnn']['stride']
        self.tstride = config['tstride']['tstride']
        self.alphabet = config['basecaller']['alphabet']
        self.qbias = config['qscore']['bias']
        self.qscale = config['qscore']['scale']
        
        #encoder
        self.gru_num_layers = config['encoder']['num_layers']
        self.gru_features = config['encoder']['features']
        self.rnn = config['encoder']['rnn']
        
        self.activation = activation
 

        # Initial convolution
        self.features = nn.Sequential(OrderedDict([
            
            ('conv0', nn.Conv1d(1, 4, kernel_size=5, stride=1, padding=5//2, bias=False)),
            ('norm0', nn.BatchNorm1d(4)),
            ('relu0',activation()),
            ('conv1',nn.Conv1d(4, 16, kernel_size=5, stride=1, padding=5//2, bias=False)),
            ('norm1', nn.BatchNorm1d(16)),
            ('relu1',activation()),
            ('conv2', nn.Conv1d(16, self.num_init_features, kernel_size=self.kernel_size, stride=1,
                                padding=self.kernel_size//2, bias=False)),
            ('norm2', nn.BatchNorm1d(self.num_init_features)),
            ('relu2', activation()),
            ('pool0', nn.MaxPool1d(kernel_size=3, stride=self.stride, padding=1))

        ]))
        
        
        
        # Each denseblock
        num_features = self.num_init_features

        self.cnn_stride = np.prod(self.tstride)*self.stride
        
        for i, num_layers in enumerate(self.block_config):
            block = _DenseBlock(num_layers=num_layers, num_input_features=num_features,
                                bn_size=self.bn_size, growth_rate=self.growth_rate, 
                                drop_rate=self.drop_rate, kernel_size=self.kernel_size, kernel_size_increment=5, activation=self.activation)
            self.features.add_module('denseblock%d' % (i + 1), block)
            num_features = num_features + num_layers * self.growth_rate
            if i != len(self.block_config) - 1:
                trans = _Transition(num_input_features=num_features, num_output_features=num_features // 2, stride = self.tstride[i], activation=self.activation)
                self.features.add_module('transition%d' % (i + 1), trans)
                num_features = num_features // 2
            

        # Final batch norm
        self.features.add_module('norm5', nn.BatchNorm1d(num_features))
        self.features.add_module('act5', activation(num_features))
        
        #self.tcn = BiTemporalConvNet(num_features, [256, 256, 256, 256, 256, 256], 3)
        
        # encoder layer
        if self.rnn:
            self.encoder = StackedBidirectionalGRU(num_features, [self.gru_features]*self.gru_num_layers, self.gru_num_layers)            
            self.encoder_output_size = self.gru_features*2
        else:
            self.encoder_output_size = num_features
        
        # decoder layer

        
        self.decoder = nn.Sequential(nn.Linear(self.encoder_output_size, 5), nn.LogSoftmax(-1))
        
        cnn_total_params = sum(p.numel() for p in self.features.parameters() if p.requires_grad)
        encoder_total_params = sum(p.numel() for p in self.encoder.parameters() if p.requires_grad) if self.rnn else 0
        decoder_total_params = sum(p.numel() for p in self.decoder.parameters() if p.requires_grad)
        # print("cnn trainable parameters: ", cnn_total_params)
        # print("encoder trainable parameters: ", encoder_total_params)
        # print("decoder trainable parameters: ", decoder_total_params)
        # print("total trainable parameters: ", cnn_total_params + encoder_total_params + decoder_total_params)
        
        for m in self.modules():
            if isinstance(m, nn.Conv1d):
                nn.init.kaiming_normal_(m.weight)
            elif isinstance(m, nn.BatchNorm1d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Linear):
                nn.init.constant_(m.bias, 0)


    def forward(self, x):
        out = self.features(x) #(N,C,T)
        #out = self.tcn(out)

        if self.rnn:
            out = out.permute(2, 0, 1) #(T, N, C)
            out = self.encoder(out) #(N, T, C)

        else:
            out = out.permute(2, 0, 1) #(T, N, C)
        #print(out.shape)
        out = self.decoder(out) #(T,N,C)

        return out# ctc需要（T,N,C）

    def use_koi(self, **kwargs):
        pass
    
    def calculate_loss(self, y, p):
        """Calculates the losses for each criterion
        
        Args:
            y (tensor): tensor with labels [batch, len]
            p (tensor): tensor with predictions [len, batch, channels]
            
        Returns:
            loss (tensor): weighted sum of losses
            losses (dict): with detached values for each loss, the weighed sum is named
                global_loss
        """
        #print(222, p.shape)
        loss = self.calculate_ctc_loss(y, p)
        # no label smooth
        losses = {'loss.global': loss.item(), 'loss.ctc': loss.item()}

        return loss, losses
    
    
    def calculate_ctc_loss(self, y, p):
        """Calculates the ctc loss
        
        Args:
            y (tensor): tensor with labels [batch, len]
            p (tensor): tensor with predictions [len, batch, channels]
            
        Returns:
            loss (tensor): weighted sum of losses
        """
        
        y_len = torch.sum(y != 0, axis = 1).to(self.device)
        p_len = torch.full((p.shape[1], ), p.shape[0]).to(self.device)
        
        loss = self.criterions["ctc"](p, y, p_len, y_len)

        return loss
    
    
    def basecall(self, p, greedy = True, *args, **kwargs):
        """Decode the predictions
         
        Args:
            p (tensor): tensor with the predictions with shape [timesteps, batch, classes]
            greedy (bool): whether to decode using a greedy approach
        Returns:
            A (list) with the decoded strings
        """
        if not isinstance(p, np.ndarray):
            p = p.cpu().numpy()

        if greedy:
            return self.decode_ctc_greedy(p, *args, **kwargs)
        else:
            return self.decode_ctc_beamsearch(p, *args, **kwargs)

    def decode_ctc_greedy(self, p, qstring = False, qscale = 1.0, qbias = 1.0, collapse_repeats = True, return_path = False, *args, **kwargs):
        """Predict the bases in a greedy approach
        Args:
            p (tensor): [len, batch, classes]
            qstring (bool): whether to return the phredq scores
            qscale (float)
            qbias (float)
        """

        decoded_predictions = list()
        #print(222, p.shape)
        for i in range(p.shape[1]):
            seq, path = viterbi_search(p[:, i, :], self.alphabet, qstring = qstring, qscale = qscale, qbias = qbias, collapse_repeats = collapse_repeats)
            if return_path:
                decoded_predictions.append((seq, path))
            else:
                decoded_predictions.append(seq)

        return decoded_predictions

    def decode_ctc_beamsearch(self, p, beam_size = 5, beam_cut_threshold = 0.1, collapse_repeats = True, *args, **kwargs):

        # print(p.shape, beam_size, beam_cut_threshold)
        # print(type(p))
        decoded_predictions = list()
        for i in range(p.shape[1]):
            seq, _ = beam_search(p[:, i, :], self.alphabet, beam_size = beam_size, beam_cut_threshold = beam_cut_threshold, collapse_repeats = collapse_repeats)
            decoded_predictions.append(seq)

        return decoded_predictions
    
    def decode(self, x, beamsize=5, threshold=1e-3, qscores=False, return_path=False):
        x = x.exp().cpu().numpy().astype(np.float32)
        if beamsize == 1 or qscores:
            seq, path  = viterbi_search(x, self.alphabet, qscores, self.qscale, self.qbias)
        else:
            seq, path = beam_search(x, self.alphabet, beamsize, threshold)
        if return_path: return seq, path
        return seq