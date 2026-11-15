# Rule archive

All four mechanisms (M1–M4) have equal standing. Each is an admissible
explanation for a structurally explicable divergence, and each is admitted only
where it predicts both the presence and the direction of the observed divergence
and its precondition holds. No mechanism is privileged over the others.

- `mechanisms.csv` — M1–M4: definition, predicted pattern and precondition
  (manuscript Table S16). The M2 precondition lists its four conditions and
  states what the rule does not require; it is written in the same form as the
  preconditions of M1, M3 and M4.
- `classification.csv` — factor-by-factor application of the mechanisms,
  including the final class (manuscript Table S17).

## Decision procedure

Implemented in `../51_b5_integration_v9.R`:

1. **Within-stream ordering** — each stream is ordered by a criterion
   appropriate to what it is licensed to estimate (population: median
   attributable fraction; individual: mutually adjusted subdistribution hazard
   ratio; genetic: counted ancestry-by-subtype cells, then multivariable
   persistence, then smallest corrected p value).
2. **Divergence recording** — a divergence is recorded where two streams place
   a factor at positions differing by two or more, or where one stream flags a
   factor as not detected while another places it.
3. **Mechanism test** — a recorded divergence is classified as structurally
   explicable only where one of M1–M4 predicts both its presence and direction
   and that mechanism's precondition holds. Otherwise the factor is classified
   as unresolved.

Labels: **C** = convergence, **SED** = structurally explicable divergence,
**UD** = unresolved divergence.
