## Chinese Remainder Theorem

A tiny, zero-dependency Python library that solves systems of simultaneous congruences with pairwise-coprime moduli, returning the unique solution modulo the product of the moduli via the extended Euclidean algorithm.

```python
from chinese_remainder_theorem import crt

# x == 2 (mod 3), x == 3 (mod 5), x == 2 (mod 7)  =>  x == 23 (mod 105)
x = crt([2, 3, 2], [3, 5, 7])
print(x)  # 23
```

Exported names: `crt(residues, moduli)` and the exception `CRTInputError` (a subclass of `ValueError`).

### Why this exists

The CRT is the kind of thing you re-derive from scratch every time you need it, usually incorrectly on the first pass. This library exists to provide one correct, dependency-free implementation you can call without thinking about Bézout coefficients or modular inverses. The trade-off: it only handles the **pairwise-coprime** case. It does not attempt to detect or resolve consistent systems whose moduli share a factor (where a solution may still exist); such inputs raise `CRTInputError`. One clear rule beats clever-but-ambiguous behaviour.

### Edge cases worth knowing

- Residues may be negative or larger than their modulus; they are normalised internally.
- A single congruence is the trivial case and returns `residue % modulus`.
- `bool` is rejected as a modulus or residue despite being a subclass of `int` — `crt([0], [True])` is almost certainly a caller mistake.
- The result is always in `[0, product_of_moduli)`.

### Running the tests

```
PYTHONPATH=src python -m unittest discover -s tests
```

## Performance

The window keeps a bounded buffer, so `push` is constant time and memory does not
grow with the length of the stream. `peak` and `trough` are linear in the window
size, which is the trade that keeps `push` cheap.

## Design notes

The window stores values eagerly rather than keeping running aggregates. Running
sums drift with floating point over long streams, and recomputing from a small
buffer is cheap enough that the drift is not worth the speed.

## Limitations

Values are coerced to floats, so very large integers lose precision. If you need
exact integer aggregates over a window, this is the wrong tool.

## Contributing

Issues and pull requests are welcome. Please keep the dependency list empty —
that constraint is the point of the project, not an oversight.

