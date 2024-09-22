## PyMOL mutation animation
The code below generates PyMOL animations for visualizing experimental evolution results.

### stop codon arising in corA in [Vial 5](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut5.part_001.vs.WT.part_001.igv.html)
```
fetch 5n77
hide everything, solvent
set_view (\
     0.012909776,    0.000930243,   -0.999916255,\
     0.104614682,   -0.994512737,    0.000425449,\
    -0.994429052,   -0.104611412,   -0.012936254,\
     0.000000000,    0.000000000, -352.192291260,\
   -49.109313965,   19.790733337,    9.610996246,\
   277.671081543,  426.713500977,  -20.000000000 )

color gray60
select stop, resi 62
select remaining, resi 0-61
select deletion, resi 63-316
select mg, resn mg
select pg0, resn pg0
hide everything, pg0

select deletionA, deletion and chain A and not resn mg
select remainingA, remaining and chain A
select stopA, stop and chain A

mset 1x180 
show cartoon, all; color gray80, all; show spheres, mg; color magenta, mg 
mview store, 0
mview store, 180
mview store, 10
mview store, 170

orient chain A
set_view (\
    -0.708065271,   -0.700900912,   -0.085915722,\
     0.551195920,   -0.624638438,    0.553181589,\
    -0.441391736,    0.344332248,    0.828618526,\
     0.000000000,    0.000000000, -249.804107666,\
   -48.848197937,   31.421995163,  -14.956176758,\
   196.947463989,  302.660766602,  -20.000000000 )
mview store, 60; mview store, 90


mdo 0: show cartoon, all; color gray80, all; show spheres, mg; color magenta, mg; hide everything, pg0
mdo 170: show cartoon, all; color gray80, all; show spheres, mg; color magenta, mg; ; hide everything, pg0

mdo 30: hide cartoon, chain B; hide everything, /5n77/I/B/MG
mdo 40: hide cartoon, chain E; hide everything, /5n77/N/E/MG
mdo 45: hide cartoon, chain D; hide everything, /5n77/L/D/MG
mdo 50: hide cartoon, chain C;
mdo 70: color red, stopA; color gray40, deletionA 
mdo 75: color gray80, stopA; color gray60, deletionA
mdo 80: color red, stopA; color gray40, deletionA
mdo 85: color gray80, stopA; color gray60, deletionA
mdo 90: color red, stopA; color gray40, deletionA
mdo 120: hide cartoon, deletionA; hide everything, mg
orient remainingA; 
mview store, 120; mview store, 130
```

### fragment deletion arising in corA in [Vial 3](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut3.part_001.vs.WT.part_001.igv.html)
```
fetch 5n77
hide everything, solvent
set_view (\
     0.012909776,    0.000930243,   -0.999916255,\
     0.104614682,   -0.994512737,    0.000425449,\
    -0.994429052,   -0.104611412,   -0.012936254,\
     0.000000000,    0.000000000, -352.192291260,\
   -49.109313965,   19.790733337,    9.610996246,\
   277.671081543,  426.713500977,  -20.000000000 )

color gray60
select stop, resi 203+240
select remaining, resi 0-202+241-999
select deletion, resi 203-240
select mg, resn mg
select pg0, resn pg0
hide everything, pg0

select deletionA, deletion and chain A and not resn mg
select remainingA, remaining and chain A
select stopA, stop and chain A

mset 1x180 
show cartoon, all; color gray80, all; show spheres, mg; color magenta, mg 
mview store, 0
mview store, 180
mview store, 10
mview store, 170

orient chain A
set_view (\
    -0.708065271,   -0.700900912,   -0.085915722,\
     0.551195920,   -0.624638438,    0.553181589,\
    -0.441391736,    0.344332248,    0.828618526,\
     0.000000000,    0.000000000, -249.804107666,\
   -48.848197937,   31.421995163,  -14.956176758,\
   196.947463989,  302.660766602,  -20.000000000 )
mview store, 60; mview store, 90


mdo 0: show cartoon, all; color gray80, all; show spheres, mg; color magenta, mg; hide everything, pg0
mdo 170: show cartoon, all; color gray80, all; show spheres, mg; color magenta, mg; ; hide everything, pg0

mdo 30: hide cartoon, chain B; hide everything, /5n77/I/B/MG
mdo 40: hide cartoon, chain E; hide everything, /5n77/N/E/MG
mdo 45: hide cartoon, chain D; hide everything, /5n77/L/D/MG
mdo 50: hide cartoon, chain C; hide everything, mg
mdo 70: color red, stopA; color red, deletionA 
mdo 75: color gray80, stopA; color gray60, deletionA
mdo 80: color red, stopA; color red, deletionA
mdo 85: color gray80, stopA; color gray60, deletionA
mdo 90: color red, stopA; color red, deletionA
mdo 120: hide cartoon, deletionA; hide everything, mg
orient deletionA; 
mview store, 120; mview store, 130
```
exporting:
```
mpng pymol_animations\corA_stop, width=1280, height=720
```
converting png sequence to gif:
```
magick convert -delay 5 -loop 0 -layers optimize .\pymol_animations\corA_stop*.png ./corA_stop.gif
```
