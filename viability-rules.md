# Reguły viability w onquadro-aligner

Program sprawdza długości trzech pętli (l1, l2, l3) oraz topologię (typ każdej pętli: propeller `p`, lateral `l`, diagonal `d`) i liczbę tetrad.

**Ważne**: wszystkie długości pętli są **obcinane (clamped)** do maksimum 4 przed porównaniem z tabelą. Czyli l=5 jest traktowane jako 4, l=6 jako 4, itd.

## Klasy topologii (wg Silva et al.)

Wzorzec typów pętli mapuje się na klasy topologii. Niektóre wzorce są dwuznaczne. Wzorce nieobserwowane eksperymentalnie (np. `dpl`) zwracają pusty zbiór klas.

| Wzorzec | Klasy topologii |
| ------- | --------------- |
| ppp | parallel |
| lll | chair |
| ldl | basket, basket2 |
| dpd | basket, basket2 |
| lpl | basket, basket2 |
| plp | basket, basket2 |
| pdp | basket, basket2 |
| pll | hybrid1 |
| pdl | hybrid2, hybrid3 |
| ldp | hybrid2, hybrid3 |
| ppl | hybrid2, hybrid3 |
| llp | hybrid2, hybrid3 |
| lpp | hybrid4 |
| dpl | *(nie istnieje)* |

## 0. Warunek wstępny

Liczba tetrad < 2 → `n/a`

## 1. Twarde reguły geometryczne

| Nr | Warunek | Decyzja |
|----|---------|---------|
| 1 | Pętla ≤1 nie-propeller (dowolna z trzech) | `not_viable` |
| 2 | Pętla diagonal długości <4 (dowolna z trzech) | `not_viable` |
| 3 | Liczba tetrad ≥4 i którakolwiek pętla ≤1 | `not_viable` |
| 4 | ≥2 pętle ≤1, topologia parallel | `viable` |
| 5 | ≥2 pętle ≤1, topologia nie-parallel | `not_viable` |
| 6 | Topologia d+pd lub d-pd, l1≥4 i l3≥4 (clamped) | `viable` |
| 7 | Topologia d+pd lub d-pd, inne długości | `not_viable` |

## 2. Nierozpoznana topologia

Jeśli klasa topologii jest pusta (np. wzorzec `dpl` nie istnieje w klasyfikacji Silva) → `unknown`

## 3. Tabela reguł (first-match-wins)

Sprawdzane w kolejności od góry do dołu. Pierwsza pasująca reguła określa decyzję.

| l1 | l2 | l3 | Klasa | Tetrady | Decyzja |
|----|----|----|-------|---------|---------|
| 1 | dowolne | 1 | parallel | ≥3 | `viable` |
| 1 | dowolne | 1 | parallel | =2 | `marginal` |
| 1 | dowolne | 1 | dowolna | dowolne | `not_viable` |
| 2 | 1 | 2 | parallel | ≥3 | `viable` |
| 2 | 1 | 2 | parallel | =2 | `marginal` |
| 2 | 1 | 2 | dowolna | dowolne | `not_viable` |
| 2 | 2 | 2 | parallel | ≥3 | `viable` |
| 2 | 2 | 2 | parallel | =2 | `marginal` |
| 2 | 2 | 2 | chair | dowolne | `marginal` |
| 2 | 2 | 2 | dowolna | dowolne | `not_viable` |
| 2 | 3 | 2 | chair | ≠3 | `viable` |
| 2 | 3 | 2 | chair | =3 | `marginal` |
| 2 | 3 | 2 | dowolna | dowolne | `not_viable` |
| 2 | ≥4 | 2 | basket | ≠3 | `viable` |
| 2 | ≥4 | 2 | basket | =3 | `marginal` |
| 2 | ≥4 | 2 | dowolna | dowolne | `not_viable` |
| 3 | 1 | 3 | parallel | ≥3 | `viable` |
| 3 | 1 | 3 | parallel | =2 | `marginal` |
| 3 | 1 | 3 | dowolna | dowolne | `not_viable` |
| 3 | 2 | 3 | hybrid3 | ≥3 | `viable` |
| 3 | 2 | 3 | hybrid3 | =2 | `marginal` |
| 3 | 2 | 3 | dowolna | dowolne | `not_viable` |
| 3 | 3 | 3 | hybrid1 | ≥3 | `viable` |
| 3 | 3 | 3 | hybrid1 | =2 | `marginal` |
| 3 | 3 | 3 | hybrid3 | dowolne | `marginal` |
| 3 | 3 | 3 | chair | dowolne | `marginal` |
| 3 | 3 | 3 | basket2 | dowolne | `marginal` |
| 3 | 3 | 3 | dowolna | dowolne | `not_viable` |
| 3 | ≥4 | 3 | basket | ≠3 | `viable` |
| 3 | ≥4 | 3 | basket | =3 | `marginal` |
| 3 | ≥4 | 3 | basket2 | dowolne | `marginal` |
| 3 | ≥4 | 3 | dowolna | dowolne | `not_viable` |
| ≥4 | 1 | ≥4 | basket | dowolne | `marginal` |
| ≥4 | 1 | ≥4 | dowolna | dowolne | `not_viable` |
| ≥4 | ≥4 | ≥4 | basket | dowolne | `marginal` |
| ≥4 | ≥4 | ≥4 | chair | dowolne | `marginal` |
| ≥4 | ≥4 | ≥4 | dowolna | dowolne | `not_viable` |

## 4. Reguły two-tetrad (l1 ≠ l3, tylko 2 tetrady)

| l1 | l2 | l3 | Klasa | Decyzja |
|----|----|----|-------|---------|
| ≥4 | ≥4 | ≥4 | basket, basket2 | `marginal` |
| ≥3 | ≥4 | ≥2 | basket, basket2 | `viable` |
| =2 | ≥4 | ≥2 | basket, basket2 | `marginal` |
| dowolne | dowolne | dowolne | basket, basket2 | `not_viable` |
| ≥4 | =2 | ≥4 | chair | `not_viable` |
| ≥3 | =2 | ≥3 | chair | `marginal` |
| dowolne | dowolne | dowolne | chair | `not_viable` |

## 5. Fallback propeller (wszystkie pętle propeller)

| l1 | l2 | l3 | Klasa | Tetrady | Decyzja |
|----|----|----|-------|---------|---------|
| ≥2 | ≥2 | ≥2 | parallel | =2 | `marginal` |
| ≥2 | ≥2 | ≥2 | parallel | ≠2 | `viable` |
| dowolne | dowolne | dowolne | parallel | dowolne | `not_viable` |

## 6. Fallback nie-propeller

| l1 | l2 | l3 | Klasa | Decyzja |
|----|----|----|-------|---------|
| ≥2 | ≥2 | ≥2 | dowolna rozpoznana | `marginal` |

## 7. Domyślne

Jeśli żadna reguła nie pasuje → `unknown`
