#Split target
```
./neighbor_target.py --clade_name Enterobacterales \
                     --neighbor_fp /home/mitsuki/afp/material/Enterobacterales/cgn/Enterobacterales_naive.neighbor \
                     --score_method naive \
                     --out_fp tmp.neighbor

./neighbor_all.py --clade_name Enterobacterales \
                     --score_method naive \
                     --out_fp tmp.neighbor

./neighbor_all.py --clade_name Enterobacterales \
                     --score_method independent \
                     --out_fp tmp.neighbor
```
