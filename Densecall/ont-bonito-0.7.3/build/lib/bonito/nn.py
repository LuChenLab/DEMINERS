"""
Bonito nn modules.
"""

import torch
from torch.nn import Module
from torch.nn.init import orthogonal_
from torch.nn.utils.fusion import fuse_conv_bn_eval


layers = {}


def register(layer):
    layer.name = layer.__name__.lower()
    layers[layer.name] = layer
    return layer


register(torch.nn.ReLU)
register(torch.nn.Tanh)


@register
class Linear(Module):
    def __init__(self, in_features, out_features, bias=True):
        super().__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.bias = bias
        self.linear = torch.nn.Linear(
            in_features=in_features, out_features=out_features, bias=bias
        )

    def forward(self, x):
        return self.linear(x)

    def to_dict(self, include_weights=False):
        res = {
            "in_features": self.in_features,
            "out_features": self.out_features,
            "bias": self.bias,
        }
        if include_weights:
            res['params'] = {
                'W': self.linear.weight,
                'b': self.linear.bias if self.bias is not None else []
            }
        return res


@register
class Swish(torch.nn.SiLU):
    pass


@register
class Clamp(Module):
    def __init__(self, min, max):
        super().__init__()
        self.min = min
        self.max = max

    def forward(self, x):
        return torch.clamp(x, min=self.min, max=self.max)

    def to_dict(self, include_weights=False):
        return {
            'min': self.min,
            'max': self.max
        }


@register
class Serial(torch.nn.Sequential):

    def __init__(self, sublayers):
        super().__init__(*sublayers)

    def forward(self, x, return_features=False):
        if return_features:
            fmaps = []
            for layer in self:
                x = layer(x)
                fmaps.append(x)
            return x, fmaps
        return super().forward(x)

    def to_dict(self, include_weights=False):
        return {
            'sublayers': [to_dict(layer, include_weights) for layer in self._modules.values()]
        }

    def __repr__(self):
        return torch.nn.ModuleList.__repr__(self)

@register
class Reverse(Module):

    def __init__(self, sublayers):
        super().__init__()
        self.layer = Serial(sublayers) if isinstance(sublayers, list) else sublayers

    def forward(self, x):
        return self.layer(x.flip(0)).flip(0)

    def to_dict(self, include_weights=False):
        if isinstance(self.layer, Serial):
            return self.layer.to_dict(include_weights)
        else:
            return {'sublayers': to_dict(self.layer, include_weights)}


@register
class BatchNorm(Module):

    def __init__(self, num_features, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True):
        super().__init__()
        self.bn = torch.nn.BatchNorm1d(num_features, eps, momentum, affine, track_running_stats)

    def forward(self, x):
        return self.bn(x)

    def to_dict(self, include_weights=False):
        res = {
            "num_features": self.bn.num_features,
            "eps": self.bn.eps,
            "momentum": self.bn.momentum,
            "affine": self.bn.affine,
            "track_running_stats": self.bn.track_running_stats
        }
        if include_weights:
            params = {}
            if res["affine"]:
                params["W"] = self.bn.weight
                params["b"] = self.bn.bias
            if res["track_running_stats"]:
                params["running_mean"] = self.bn.running_mean
                params["running_var"] = self.bn.running_var
            res["params"] = params
        return res


@register
class Convolution(Module):

    def __init__(self, insize, size, winlen, stride=1, padding=0, bias=True, activation=None, norm=None):
        super().__init__()
        self.conv = torch.nn.Conv1d(insize, size, winlen, stride=stride, padding=padding, bias=bias)
        self.activation = layers.get(activation, lambda: activation)()
        if isinstance(norm, dict):
            self.norm = from_dict(norm)
        elif isinstance(norm, str):
            self.norm = layers[norm](size)
        else:
            self.norm = norm

    def forward(self, x):
        h = self.conv(x)
        if self.norm is not None:
            h = self.norm(h)
        if self.activation is not None:
            h = self.activation(h)
        return h

    def to_dict(self, include_weights=False):
        res = {
            "insize": self.conv.in_channels,
            "size": self.conv.out_channels,
            "bias": self.conv.bias is not None,
            "winlen": self.conv.kernel_size[0],
            "stride": self.conv.stride[0],
            "padding": self.conv.padding[0],
        }
        if self.activation is not None:
            res["activation"] = self.activation.name
        if self.norm is not None:
            res["norm"] = to_dict(self.norm, include_weights)
            #simplify default case e.g. norm="batchnorm"
            if not include_weights and self.norm.name in layers:
                if res["norm"] == to_dict(layers[self.norm.name](res["size"])):
                    res["norm"] = self.norm.name

        if include_weights:
            res['params'] = {
                'W': self.conv.weight, 'b': self.conv.bias if self.conv.bias is not None else []
            }
        return res


@register
class LinearCRFEncoder(Module):

    def __init__(self, insize, n_base, state_len, bias=True, scale=None, activation=None, blank_score=None, expand_blanks=True, permute=None):
        super().__init__()
        self.scale = scale
        self.n_base = n_base
        self.state_len = state_len
        self.blank_score = blank_score
        self.expand_blanks = expand_blanks
        size = (n_base + 1) * n_base**state_len if blank_score is None else n_base**(state_len + 1)
        self.linear = torch.nn.Linear(insize, size, bias=bias)
        self.activation = layers.get(activation, lambda: activation)()
        self.permute = permute

    def forward(self, x):
        if self.permute is not None:
            x = x.permute(*self.permute)
        scores = self.linear(x)
        if self.activation is not None:
            scores = self.activation(scores)
        if self.scale is not None:
            scores = scores * self.scale
        if self.blank_score is not None and self.expand_blanks:
            T, N, C = scores.shape
            scores = torch.nn.functional.pad(
                scores.view(T, N, C // self.n_base, self.n_base),
                (1, 0, 0, 0, 0, 0, 0, 0),
                value=self.blank_score
            ).view(T, N, -1)
        return scores

    def to_dict(self, include_weights=False):
        res = {
            'insize': self.linear.in_features,
            'n_base': self.n_base,
            'state_len': self.state_len,
            'bias': self.linear.bias is not None,
            'scale': self.scale,
            'blank_score': self.blank_score,
            'expand_blanks': self.expand_blanks,
        }
        if self.activation is not None:
            res['activation'] = self.activation.name
        if self.permute is not None:
            res['permute'] = self.permute
        if include_weights:
            res['params'] = {
                'W': self.linear.weight, 'b': self.linear.bias
                if self.linear.bias is not None else []
            }
        return res

    def extra_repr(self):
        rep = 'n_base={}, state_len={}, scale={}, blank_score={}, expand_blanks={}'.format(
            self.n_base, self.state_len, self.scale, self.blank_score, self.expand_blanks
        )
        if self.permute:
            rep += f', permute={self.permute}'
        return rep
            

@register
class Permute(Module):

    def __init__(self, dims):
        super().__init__()
        self.dims = dims

    def forward(self, x):
        return x.permute(*self.dims)

    def to_dict(self, include_weights=False):
        return {'dims': self.dims}

    def extra_repr(self):
        return 'dims={}'.format(self.dims)


def truncated_normal(size, dtype=torch.float32, device=None, num_resample=5):
    x = torch.empty(size + (num_resample,), dtype=torch.float32, device=device).normal_()
    i = ((x < 2) & (x > -2)).max(-1, keepdim=True)[1]
    return torch.clamp_(x.gather(-1, i).squeeze(-1), -2, 2)


class RNNWrapper(Module):
    def __init__(
            self, rnn_type, *args, reverse=False, orthogonal_weight_init=True, disable_state_bias=True, bidirectional=False, **kwargs
    ):
        super().__init__()
        if reverse and bidirectional:
            raise Exception("'reverse' and 'bidirectional' should not both be set to True")
        self.reverse = reverse
        self.rnn = rnn_type(*args, bidirectional=bidirectional, **kwargs)
        self.init_orthogonal(orthogonal_weight_init)
        self.init_biases()
        if disable_state_bias: self.disable_state_bias()

    def forward(self, x):
        if self.reverse: x = x.flip(0)
        y, h = self.rnn(x)
        if self.reverse: y = y.flip(0)
        return y

    def init_biases(self, types=('bias_ih',)):
        for name, param in self.rnn.named_parameters():
            if any(k in name for k in types):
                with torch.no_grad():
                    param.set_(0.5*truncated_normal(param.shape, dtype=param.dtype, device=param.device))

    def init_orthogonal(self, types=True):
        if not types: return
        if types == True: types = ('weight_ih', 'weight_hh')
        for name, x in self.rnn.named_parameters():
            if any(k in name for k in types):
                for i in range(0, x.size(0), self.rnn.hidden_size):
                    orthogonal_(x[i:i+self.rnn.hidden_size])

    def disable_state_bias(self):
        for name, x in self.rnn.named_parameters():
            if 'bias_hh' in name:
                x.requires_grad = False
                x.zero_()

    def extra_repr(self):
        return 'reverse={}'.format(bool(self.reverse))


@register
class LSTM(RNNWrapper):

    def __init__(self, size, insize, bias=True, reverse=False):
        super().__init__(torch.nn.LSTM, insize, size, bias=bias, reverse=reverse)

    def to_dict(self, include_weights=False):
        res = {
            'size': self.rnn.hidden_size,
            'insize': self.rnn.input_size,
            'bias': self.rnn.bias,
            'reverse': self.reverse,
        }
        if include_weights:
            res['params'] = {
                'iW': self.rnn.weight_ih_l0.reshape(4, self.rnn.hidden_size, self.rnn.input_size),
                'sW': self.rnn.weight_hh_l0.reshape(4, self.rnn.hidden_size, self.rnn.hidden_size),
                'b': self.rnn.bias_ih_l0.reshape(4, self.rnn.hidden_size)
            }
        return res


def to_dict(layer, include_weights=False):
    if hasattr(layer, 'to_dict'):
        return {'type': layer.name, **layer.to_dict(include_weights)}
    return {'type': layer.name}


def from_dict(model_dict, layer_types=None):
    if not isinstance(model_dict, dict):
        # enable model_dict to contain concrete objects, handy in nb creation etc
        return model_dict
    model_dict = model_dict.copy()
    if layer_types is None:
        layer_types = layers
    type_name = model_dict.pop('type')
    typ = layer_types[type_name]
    if hasattr(typ, "from_dict"):
        return typ.from_dict(model_dict, layer_types)
    if 'sublayers' in model_dict:
        sublayers = model_dict['sublayers']
        model_dict['sublayers'] = [
            from_dict(x, layer_types) for x in sublayers
        ] if isinstance(sublayers, list) else from_dict(sublayers, layer_types)
    try:
        layer = typ(**model_dict)
    except Exception as e:
        raise Exception(f'Failed to build layer of type {typ} with args {model_dict}') from e
    return layer


def fuse_bn_(m):
    """
    Sets the module m to eval mode and if a Convolution fuses any batchnorm layer.
    """
    m.training = False
    if isinstance(m, Convolution) and isinstance(m.norm, BatchNorm):
        m.conv = fuse_conv_bn_eval(m.conv, m.norm.bn)
        m.norm = None
