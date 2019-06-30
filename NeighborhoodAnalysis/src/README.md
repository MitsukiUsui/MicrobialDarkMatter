```
./neighbor_target.py --clade_name Enterobacterales \
                     --score_method naive \
                     --neighbor_fp /home/mitsuki/afp/material/Enterobacterales/cgn/Enterobacterales_naive.neighbor \
                     --out_fp tmp.neighbor

./neighbor_all.py --clade_name Enterobacterales \
                     --score_method naive \
                     --out_fp tmp.neighbor

./neighbor_all.py --clade_name Enterobacterales \
                     --score_method independent \
                     --out_fp tmp.neighbor

./neighbor_all.py --clade_name Enterobacterales \
                     --score_method conditional \
                     --out_fp tmp.neighbor
```
