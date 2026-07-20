Thanks for taking the time to contribute — we appreciate it very much!

We welcome contributions to the core engine of ReactiveDynamics.jl — refining and extending the modeling metalanguage, the resource allocator, and the analysis/observability layer for improved user experience and modeling expressivity — as well as to the documentation and the worked decision case studies.

Please start with the design records under [`spec/`](spec), which document the engine's behavior and the reasoning behind it:

- [`spec/STATUS.md`](spec/STATUS.md) — the single "what is the state, what is left" index. Start here.
- [`spec/CONTRACT_DRAFT.md`](spec/CONTRACT_DRAFT.md) — the normative operational-semantics specification (§1–§15).
- [`spec/adr/`](spec/adr) — the Architecture Decision Records: *why* the engine is the way it is.
- [`spec/INVENTORY.md`](spec/INVENTORY.md) — the current-source map (module map, public-API audit, stepping trace, AlgebraicAgents touchpoints).

Two conventions matter. Decisions are changed by **adding a new ADR**, not by rewriting an existing one — the record is append-only. And when documentation and code disagree, **the code is the source of truth**; fix the docs to match, or open an issue if the code looks wrong.

Before submitting a pull request, please run the semantic test suite and the formatter locally (there is no CI in this repo):

```bash
julia --project=. -e 'using Pkg; Pkg.test()'          # the semantic suite
julia -m Runic --check .                               # formatting (Runic — opinionated, no config)
```

[Runic](https://github.com/fredrikekre/Runic.jl) is non-configurable; `julia -m Runic --inplace src test ext dev docs/make.jl demo` applies the formatting in place.

Contributions to ReactiveDynamics.jl are welcome in the following forms:

- Modifying the code or documentation via a pull request.
- Reporting bugs or suggesting enhancements in the project's [GitHub Issues](https://github.com/Merck/ReactiveDynamics.jl/issues). If you propose a feature for future development, we are happy to discuss and take on the implementation.
