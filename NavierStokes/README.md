# NavierStokes Command Line
## Mode
Alias|Mode
---|---
-td|TDMA solver
### Optional Variables
#### TDMA solver
Variable|Meaning
---|---
-aa|diagonal column below the middle
-bb|middle diagonal column value
-cc|diagonal column above the middle
-dd|the right column (b)

#### Example
```bash
NavierStokes -aa ARRAY[-5,-5,-5,-5] -bb ARRAY[20,15,15,15,10] -cc ARRAY[-5,-5,-5,-5] -dd ARRAY[1100,100,100,100] -td 
```
