# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import argparse
import pathlib
import string
import os
import torch
from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained, MSATransformer
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import itertools
from typing import List, Tuple
import numpy as np

def remove_insertions(sequence: str) -> str:
    deletekeys = dict.fromkeys(string.ascii_lowercase)
    deletekeys["."] = None
    deletekeys["*"] = None
    translation = str.maketrans(deletekeys)
    return sequence.translate(translation)

def read_msa(filename: str, nseq: int) -> List[Tuple[str, str]]:
    msa = [
        (record.description, remove_insertions(str(record.seq)))
        for record in itertools.islice(SeqIO.parse(filename, "fasta"), nseq)
    ]
    return msa

def create_parser():
    parser = argparse.ArgumentParser(description="Batch mutation scoring with ESM-1v models.")
    parser.add_argument("--model-location", type=str, nargs="+", help="ESM model paths")
    parser.add_argument("--dms-input", type=pathlib.Path, help="CSV file with sequences and mutations")
    parser.add_argument("--mutation-col", type=str, default="mutant", help="Mutation column")
    parser.add_argument("--dms-output", type=pathlib.Path, help="CSV output file")
    parser.add_argument("--offset-idx", type=int, default=0, help="Offset to subtract from mutation index")
    parser.add_argument("--scoring-strategy", type=str, default="wt-marginals", choices=["wt-marginals", "pseudo-ppl", "masked-marginals"])
    parser.add_argument("--msa-path", type=pathlib.Path, help="MSA file path if using MSA Transformer")
    parser.add_argument("--msa-samples", type=int, default=400)
    parser.add_argument("--nogpu", action="store_true")
    return parser

def label_row(row, token_probs, alphabet):
    wt, idx, mt = row["mutant"][0], int(row["mutant"][1:-1]) - int(row["offset"]), row["mutant"][-1]
    seq = row["sequence"]
    assert seq[idx] == wt, f"Wildtype mismatch at {idx}: expected {wt}, got {seq[idx]}"
    wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)
    score = token_probs[0, 1 + idx, mt_encoded] - token_probs[0, 1 + idx, wt_encoded]
    return score.item()

def compute_pppl(row, model, alphabet):
    wt, idx, mt = row["mutant"][0], int(row["mutant"][1:-1]) - int(row["offset"]), row["mutant"][-1]
    seq = row["sequence"]
    assert seq[idx] == wt, f"Wildtype mismatch at {idx}: expected {wt}, got {seq[idx]}"
    sequence = seq[:idx] + mt + seq[idx + 1:]
    data = [("protein", sequence)]
    batch_converter = alphabet.get_batch_converter()
    _, _, batch_tokens = batch_converter(data)
    log_probs = []
    for i in range(1, len(sequence) - 1):
        batch_tokens_masked = batch_tokens.clone()
        batch_tokens_masked[0, i] = alphabet.mask_idx
        with torch.no_grad():
            token_probs = torch.log_softmax(model(batch_tokens_masked.cuda())['logits'], dim=-1)
        log_probs.append(token_probs[0, i, alphabet.get_idx(sequence[i])].item())
    return sum(log_probs)

def select_least_loaded_gpu():
    try:
        import pynvml
        pynvml.nvmlInit()
        min_used = float("inf")
        best_gpu = 0
        for i in range(torch.cuda.device_count()):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            mem = pynvml.nvmlDeviceGetMemoryInfo(handle)
            used = mem.used / mem.total
            if used < min_used:
                min_used = used
                best_gpu = i
        pynvml.nvmlShutdown()
        return best_gpu
    except Exception as e:
        print("Could not auto-select GPU:", e)
        return 0

def main(args):
    df = pd.read_csv(args.dms_input)
    df.columns = df.columns.str.strip()
    df = df.rename(columns={"sequence": "sequence", "mutant": "mutant", "offset": "offset"})

    for model_location in args.model_location:
        model, alphabet = pretrained.load_model_and_alphabet(model_location)
        model.eval()

        if torch.cuda.is_available() and not args.nogpu:
            gpu_idx = select_least_loaded_gpu()
            torch.cuda.set_device(gpu_idx)
            model = model.cuda()
            print(f"Transferred model to GPU:{gpu_idx} -> {model_location}")

        batch_converter = alphabet.get_batch_converter()

        if isinstance(model, MSATransformer):
            data = [read_msa(args.msa_path, args.msa_samples)]
            assert args.scoring_strategy == "masked-marginals"
            _, _, batch_tokens = batch_converter(data)
            all_token_probs = []
            for i in tqdm(range(batch_tokens.size(2))):
                batch_tokens_masked = batch_tokens.clone()
                batch_tokens_masked[0, 0, i] = alphabet.mask_idx
                with torch.no_grad():
                    token_probs = torch.log_softmax(model(batch_tokens_masked.cuda())['logits'], dim=-1)
                all_token_probs.append(token_probs[:, 0, i])
            token_probs = torch.cat(all_token_probs, dim=0).unsqueeze(0)
            df[model_location] = df.apply(lambda row: label_row(row, token_probs, alphabet), axis=1)

        else:
            for seq in tqdm(df["sequence"].unique(), desc=f"Scoring {model_location}"):
                subset = df[df["sequence"] == seq]
                data = [("protein", seq)]
                _, _, batch_tokens = batch_converter(data)
                with torch.no_grad():
                    token_probs = torch.log_softmax(model(batch_tokens.cuda())['logits'], dim=-1)
                df.loc[subset.index, model_location] = subset.apply(
                    lambda row: label_row(row, token_probs, alphabet), axis=1
                )

    df.to_csv(args.dms_output, index=False)

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)
