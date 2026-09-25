# Adaptation of the requested SLM scaffold

Source: [ChaitanyaK77/Building-a-Small-Language-Model-SLM-](https://github.com/ChaitanyaK77/Building-a-Small-Language-Model-SLM-), inspected at commit `243021c36824aba864c32e15b591044318d2216f` (2025-06-07). The reference notebook is `Small_Language_Model_Building.ipynb`. The upstream MIT notice is retained in `references/SLM-LICENSE`.

This project follows its setup sequence and adapts its small decoder architecture. It is organized as importable modules and command-line stages so that runs are repeatable without depending on notebook execution order.

| Scaffold component | Project adaptation |
|---|---|
| TinyStories loading | IMPI-guided, provenance-tracked retrieval of real human CDS |
| GPT-2 BPE tokenizer | Deterministic four-base and exhaustive 64-codon encoders |
| Binary token arrays | uint8 memory-mapped arrays with explicit per-sequence offsets |
| Random windows in concatenated text | Whole-CDS batches; no biological boundary crossing |
| Learned position table | Fixed sinusoidal positions measured in nucleotide coordinates |
| Causal attention, pre-norm blocks, GELU MLP | Retained in a compact configurable decoder |
| Weight-tied embedding/output head | Retained; vocabulary sizes are derived from the tokenizer |
| AdamW, warmup, cosine decay, clipping | Validated schedule updated once per optimizer update |
| Mixed precision | float32 default; optional supported CUDA bfloat16 |
| Best validation checkpoint | Final/scheduled resumable checkpoints; validation explicitly deferred |
| Generated story inspection | Training diagnostics and simple empirical baselines only |

Several notebook details are intentionally corrected for this application: the learning-rate floor cannot exceed the peak; the schedule is tied to optimizer updates rather than accumulation microsteps; gradients are normalized by valid target counts; sequence storage cannot create transitions across genes; seeds and optimizer/RNG state are saved; masking does not enlarge the requested vocabulary. The notebook's large text-model dimensions and GPU-memory reservation are not copied into this small biological dataset.

The code is not a claim that language-model performance transfers to protein function. The scaffold supplies an engineering structure; the biological sampling and interpretation require their own protocol.
