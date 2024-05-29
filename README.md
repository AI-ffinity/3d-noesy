# Idea

Theoretically, the closer the two nuclei together, (here: ¹H) the stronger the NOE signal between them.  
It it practically true though?

## Hypotheses to test
- **H1**: the most intense peak per spin system is $H^A_{i}$
- **H2**: the 2nd most intense is $H^A_{i-1}$ and the 3rd is $H^N_{i-1}$
- **H3**: the most intense inter-residual $H^N$ peak per spin system is $H^N_{i-1}$

---

# Results
> None of them holds true

### Examples

**2K52 (74 residues)**

| Atom name   |   1st highest |   2nd highest |   3rd highest |   4th or lower |
|:------------|--------------:|--------------:|--------------:|---------------:|
| HA_i        |             7 |            32 |            17 |             18 |
| HA_i-1      |            38 |            13 |             6 |             17 |
| H_i-1       |             9 |             4 |            16 |             35 |

**2LEA (100 residues)**

| Atom name   |   1st highest |   2nd highest |   3rd highest |   4th or lower |
|:------------|--------------:|--------------:|--------------:|---------------:|
| HA_i        |            30 |            38 |            20 |              5 |
| HA_i-1      |            41 |            14 |            10 |             25 |
| H_i-1       |             7 |            15 |            18 |             43 |

**2LTM (107 residues)**

| Atom name   |   1st highest |   2nd highest |   3rd highest |   4th or lower |
|:------------|--------------:|--------------:|--------------:|---------------:|
| HA_i        |            14 |            40 |            22 |             10 |
| HA_i-1      |            31 |            13 |            10 |             31 |
| H_i-1       |            18 |             9 |             9 |             33 |
