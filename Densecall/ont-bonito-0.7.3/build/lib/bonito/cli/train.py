#!/usr/bin/env python3

"""
Bonito training.
"""

import os
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from pathlib import Path
from importlib import import_module

from bonito.data import load_numpy, load_script
from bonito.util import __models__, default_config, default_data
from bonito.util import load_model, load_symbol, init, half_supported
from bonito.training import load_state, Trainer

import toml
import torch
import numpy as np
from torch.utils.data import DataLoader


def main(args):

    workdir = os.path.expanduser(args.training_directory)

    if os.path.exists(workdir) and not args.force:
        print("[error] %s exists, use -f to force continue training." % workdir)
        exit(1)

    init(args.seed, args.device, (not args.nondeterministic))
    device = torch.device(args.device)

    if not args.pretrained:
        config = toml.load(args.config)
    else:
        dirname = args.pretrained
        if not os.path.isdir(dirname) and os.path.isdir(os.path.join(__models__, dirname)):
            dirname = os.path.join(__models__, dirname)
        pretrain_file = os.path.join(dirname, 'config.toml')
        config = toml.load(pretrain_file)
        if 'lr_scheduler' in config:
            print(f"[ignoring 'lr_scheduler' in --pretrained config]")
            del config['lr_scheduler']

    argsdict = dict(training=vars(args))

    print("[loading model]")
    if args.pretrained:
        print("[using pretrained model {}]".format(args.pretrained))
        model = load_model(args.pretrained, device, half=False)
    else:
        model = load_symbol(config, 'Model')(config)

    print("[loading data]")
    try:
        train_loader_kwargs, valid_loader_kwargs = load_numpy(
            args.chunks, args.directory, valid_chunks = args.valid_chunks
        )
    except FileNotFoundError:
        train_loader_kwargs, valid_loader_kwargs = load_script(
            args.directory,
            seed=args.seed,
            chunks=args.chunks,
            valid_chunks=args.valid_chunks,
            n_pre_context_bases=getattr(model, "n_pre_context_bases", 0),
            n_post_context_bases=getattr(model, "n_post_context_bases", 0),
        )

    loader_kwargs = {
        "batch_size": args.batch, "num_workers": 4, "pin_memory": True
    }
    train_loader = DataLoader(**loader_kwargs, **train_loader_kwargs)
    valid_loader = DataLoader(**loader_kwargs, **valid_loader_kwargs)

    os.makedirs(workdir, exist_ok=True)
    toml.dump({**config, **argsdict}, open(os.path.join(workdir, 'config.toml'), 'w'))

    if config.get("lr_scheduler"):
        sched_config = config["lr_scheduler"]
        lr_scheduler_fn = getattr(
            import_module(sched_config["package"]), sched_config["symbol"]
        )(**sched_config)
    else:
        lr_scheduler_fn = None

    trainer = Trainer(
        model, device, train_loader, valid_loader,
        use_amp=half_supported() and not args.no_amp,
        lr_scheduler_fn=lr_scheduler_fn,
        restore_optim=args.restore_optim,
        save_optim_every=args.save_optim_every,
        grad_accum_split=args.grad_accum_split,
        quantile_grad_clip=args.quantile_grad_clip
    )

    if (',' in args.lr):
        lr = [float(x) for x in args.lr.split(',')]
    else:
        lr = float(args.lr)
    trainer.fit(workdir, args.epochs, lr)

def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument("training_directory")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--config', default=default_config)
    group.add_argument('--pretrained', default="")
    parser.add_argument("--directory", type=Path)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--lr", default='2e-3')
    parser.add_argument("--seed", default=25, type=int)
    parser.add_argument("--epochs", default=5, type=int)
    parser.add_argument("--batch", default=64, type=int)
    parser.add_argument("--chunks", default=0, type=int)
    parser.add_argument("--valid-chunks", default=None, type=int)
    parser.add_argument("--no-amp", action="store_true", default=False)
    parser.add_argument("-f", "--force", action="store_true", default=False)
    parser.add_argument("--restore-optim", action="store_true", default=False)
    parser.add_argument("--nondeterministic", action="store_true", default=False)
    parser.add_argument("--save-optim-every", default=10, type=int)
    parser.add_argument("--grad-accum-split", default=1, type=int)
    parser.add_argument("--quantile-grad-clip", action="store_true", default=False)
    return parser
