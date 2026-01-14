# Get data frame of transition states from Biogeographical Stochastic Mapping

Get data frame of transition states from Biogeographical Stochastic
Mapping

## Usage

``` r
get_insert_df(bsm, phyllip.file, max.range.size, include_null_range = TRUE)
```

## Arguments

- bsm:

  results from the function
  [`calc_bsm()`](https://gabrielnakamura.github.io/Herodotools/reference/calc_bsm.md)

- phyllip.file:

  path to the phyllip file used in the original model

- max.range.size:

  maximum range size used in the biogeographical reconstruction

- include_null_range:

  default TRUE.

## Value

a list with same length as the number of stochastic mapping. Each list
element is a data frame with columns:

- `child` = child node number,

- `parent` = ancestor node number,

- `event_time` = time of the event from parent to child,

- `node_area` = the range area in the event,

- `event_type` = the type of event,

- `source` = can be 'cladogenesis' or 'anagenesis'
