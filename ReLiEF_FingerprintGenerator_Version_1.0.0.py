##############################################
# This script creates Receptor Ligand surfacE Feature (ReLiEF) fingerprints from a directory of *.wrl files corresponding to the sixty conformations belonging to Abl protein kinase and probed by nuclear magnetic resonance
# as detailed in the following publication:
# Xie, T.; Saleh, T.; Rossi, P.; Kalodimos, C. G.; Conformational states dynamically populated by a kinase determine its function, Science  2020, 370 (6513).  DOI: 10.1126/science.abc2754
# 
# Copyright 2023 Benjamin M. Samudio
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Ben Samudio, May 2023
# Towards Alleviating Suffering
###############################################




import re
import math
import statistics
import os


distance_codes = [[0,0.1,'A'],[0.1,0.2,'B'],[0.2,0.3,'C'],[0.3,0.4,'D'],[0.4,0.5,'E'],[0.5,0.6,'F'],[0.6,0.7,'G'],[0.7,0.8,'H'],[0.8,0.9,'I'],[0.9,1,'J'],[1,1.1,'K'],[1.1,1.2,'L'],[1.2,1.3,'M'],[1.3,1.4,'N'],[1.4,1.5,'O'],[1.5,1.6,'P'],[1.6,1.7,'Q'],[1.7,1.8,'R'],[1.8,1.9,'S'],[1.9,2,'T'],[2,2.1,'U'],[2.1,2.2,'V'],[2.2,2.3,'W'],[2.3,2.4,'X'],[2.4,2.5,'Y'],[2.5,2.6,'Z'],[2.6,2.7,'AA'],[2.7,2.8,'AB'],[2.8,2.9,'AC'],[2.9,3,'AD'],[3,3.1,'AE'],[3.1,3.2,'AF'],[3.2,3.3,'AG'],[3.3,3.4,'AH'],[3.4,3.5,'AI'],[3.5,3.6,'AJ'],[3.6,3.7,'AK'],[3.7,3.8,'AL'],[3.8,3.9,'AM'],[3.9,4,'AN'],[4,4.1,'AO'],[4.1,4.2,'AP'],[4.2,4.3,'AQ'],[4.3,4.4,'AR'],[4.4,4.5,'AS'],[4.5,4.6,'AT'],[4.6,4.7,'AU'],[4.7,4.8,'AV'],[4.8,4.9,'AW'],[4.9,5,'AX'],[5,5.1,'AY'],[5.1,5.2,'AZ'],[5.2,5.3,'BA'],[5.3,5.4,'BB'],[5.4,5.5,'BC'],[5.5,5.6,'BD'],[5.6,5.7,'BE'],[5.7,5.8,'BF'],[5.8,5.9,'BG'],[5.9,6,'BH'],[6,6.1,'BI'],[6.1,6.2,'BJ'],[6.2,6.3,'BK'],[6.3,6.4,'BL'],[6.4,6.5,'BM'],[6.5,6.6,'BN'],[6.6,6.7,'BO'],[6.7,6.8,'BP'],[6.8,6.9,'BQ'],[6.9,7,'BR'],[7,7.1,'BS'],[7.1,7.2,'BT'],[7.2,7.3,'BU'],[7.3,7.4,'BV'],[7.4,7.5,'BW'],[7.5,7.6,'BX'],[7.6,7.7,'BY'],[7.7,7.8,'BZ'],[7.8,7.9,'CA'],[7.9,8,'CB'],[8,8.1,'CC'],[8.1,8.2,'CD'],[8.2,8.3,'CE'],[8.3,8.4,'CF'],[8.4,8.5,'CG'],[8.5,8.6,'CH'],[8.6,8.7,'CI'],[8.7,8.8,'CJ'],[8.8,8.9,'CK'],[8.9,9,'CL'],[9,9.1,'CM'],[9.1,9.2,'CN'],[9.2,9.3,'CO'],[9.3,9.4,'CP'],[9.4,9.5,'CQ'],[9.5,9.6,'CR'],[9.6,9.7,'CS'],[9.7,9.8,'CT'],[9.8,9.9,'CU'],[9.9,10,'CV'],[10,10.1,'CW'],[10.1,10.2,'CX'],[10.2,10.3,'CY'],[10.3,10.4,'CZ'],[10.4,10.5,'DA'],[10.5,10.6,'DB'],[10.6,10.7,'DC'],[10.7,10.8,'DD'],[10.8,10.9,'DE'],[10.9,11,'DF'],[11,11.1,'DG'],[11.1,11.2,'DH'],[11.2,11.3,'DI'],[11.3,11.4,'DJ'],[11.4,11.5,'DK'],[11.5,11.6,'DL'],[11.6,11.7,'DM'],[11.7,11.8,'DN'],[11.8,11.9,'DO'],[11.9,12,'DP'],[12,12.1,'DQ'],[12.1,12.2,'DR'],[12.2,12.3,'DS'],[12.3,12.4,'DT'],[12.4,12.5,'DU'],[12.5,12.6,'DV'],[12.6,12.7,'DW'],[12.7,12.8,'DX'],[12.8,12.9,'DY'],[12.9,13,'DZ'],[13,13.1,'EA'],[13.1,13.2,'EB'],[13.2,13.3,'EC'],[13.3,13.4,'ED'],[13.4,13.5,'EE'],[13.5,13.6,'EF'],[13.6,13.7,'EG'],[13.7,13.8,'EH'],[13.8,13.9,'EI'],[13.9,14,'EJ'],[14,14.1,'EK'],[14.1,14.2,'EL'],[14.2,14.3,'EM'],[14.3,14.4,'EN'],[14.4,14.5,'EO'],[14.5,14.6,'EP'],[14.6,14.7,'EQ'],[14.7,14.8,'ER'],[14.8,14.9,'ES'],[14.9,15,'ET'],[15,15.1,'EU'],[15.1,15.2,'EV'],[15.2,15.3,'EW'],[15.3,15.4,'EX'],[15.4,15.5,'EY'],[15.5,15.6,'EZ'],[15.6,15.7,'FA'],[15.7,15.8,'FB'],[15.8,15.9,'FC'],[15.9,16,'FD'],[16,16.1,'FE'],[16.1,16.2,'FF'],[16.2,16.3,'FG'],[16.3,16.4,'FH'],[16.4,16.5,'FI'],[16.5,16.6,'FJ'],[16.6,16.7,'FK'],[16.7,16.8,'FL'],[16.8,16.9,'FM'],[16.9,17,'FN'],[17,17.1,'FO'],[17.1,17.2,'FP'],[17.2,17.3,'FQ'],[17.3,17.4,'FR'],[17.4,17.5,'FS'],[17.5,17.6,'FT'],[17.6,17.7,'FU'],[17.7,17.8,'FV'],[17.8,17.9,'FW'],[17.9,18,'FX'],[18,18.1,'FY'],[18.1,18.2,'FZ'],[18.2,18.3,'GA'],[18.3,18.4,'GB'],[18.4,18.5,'GC'],[18.5,18.6,'GD'],[18.6,18.7,'GE'],[18.7,18.8,'GF'],[18.8,18.9,'GG'],[18.9,19,'GH'],[19,19.1,'GI'],[19.1,19.2,'GJ'],[19.2,19.3,'GK'],[19.3,19.4,'GL'],[19.4,19.5,'GM'],[19.5,19.6,'GN'],[19.6,19.7,'GO'],[19.7,19.8,'GP'],[19.8,19.9,'GQ'],[19.9,20,'GR'],[20,20.1,'GS'],[20.1,20.2,'GT'],[20.2,20.3,'GU'],[20.3,20.4,'GV'],[20.4,20.5,'GW'],[20.5,20.6,'GX'],[20.6,20.7,'GY'],[20.7,20.8,'GZ'],[20.8,20.9,'HA'],[20.9,21,'HB'],[21,21.1,'HC'],[21.1,21.2,'HD'],[21.2,21.3,'HE'],[21.3,21.4,'HF'],[21.4,21.5,'HG'],[21.5,21.6,'HH'],[21.6,21.7,'HI'],[21.7,21.8,'HJ'],[21.8,21.9,'HK'],[21.9,22,'HL'],[22,22.1,'HM'],[22.1,22.2,'HN'],[22.2,22.3,'HO'],[22.3,22.4,'HP'],[22.4,22.5,'HQ'],[22.5,22.6,'HR'],[22.6,22.7,'HS'],[22.7,22.8,'HT'],[22.8,22.9,'HU'],[22.9,23,'HV'],[23,23.1,'HW'],[23.1,23.2,'HX'],[23.2,23.3,'HY'],[23.3,23.4,'HZ'],[23.4,23.5,'IA'],[23.5,23.6,'IB'],[23.6,23.7,'IC'],[23.7,23.8,'ID'],[23.8,23.9,'IE'],[23.9,24,'IF'],[24,24.1,'IG'],[24.1,24.2,'IH'],[24.2,24.3,'II'],[24.3,24.4,'IJ'],[24.4,24.5,'IK'],[24.5,24.6,'IL'],[24.6,24.7,'IM'],[24.7,24.8,'IN'],[24.8,24.9,'IO'],[24.9,25,'IP'],[25,25.1,'IQ'],[25.1,25.2,'IR'],[25.2,25.3,'IS'],[25.3,25.4,'IT'],[25.4,25.5,'IU'],[25.5,25.6,'IV'],[25.6,25.7,'IW'],[25.7,25.8,'IX'],[25.8,25.9,'IY'],[25.9,26,'IZ'],[26,26.1,'JA'],[26.1,26.2,'JB'],[26.2,26.3,'JC'],[26.3,26.4,'JD'],[26.4,26.5,'JE'],[26.5,26.6,'JF'],[26.6,26.7,'JG'],[26.7,26.8,'JH'],[26.8,26.9,'JI'],[26.9,27,'JJ'],[27,27.1,'JK'],[27.1,27.2,'JL'],[27.2,27.3,'JM'],[27.3,27.4,'JN'],[27.4,27.5,'JO'],[27.5,27.6,'JP'],[27.6,27.7,'JQ'],[27.7,27.8,'JR'],[27.8,27.9,'JS'],[27.9,28,'JT'],[28,28.1,'JU'],[28.1,28.2,'JV'],[28.2,28.3,'JW'],[28.3,28.4,'JX'],[28.4,28.5,'JY'],[28.5,28.6,'JZ'],[28.6,28.7,'KA'],[28.7,28.8,'KB'],[28.8,28.9,'KC'],[28.9,29,'KD'],[29,29.1,'KE'],[29.1,29.2,'KF'],[29.2,29.3,'KG'],[29.3,29.4,'KH'],[29.4,29.5,'KI'],[29.5,29.6,'KJ'],[29.6,29.7,'KK'],[29.7,29.8,'KL'],[29.8,29.9,'KM'],[29.9,30,'KN'],[30,30.1,'KO'],[30.1,30.2,'KP'],[30.2,30.3,'KQ'],[30.3,30.4,'KR'],[30.4,30.5,'KS'],[30.5,30.6,'KT'],[30.6,30.7,'KU'],[30.7,30.8,'KV'],[30.8,30.9,'KW'],[30.9,31,'KX'],[31,31.1,'KY'],[31.1,31.2,'KZ'],[31.2,31.3,'LA'],[31.3,31.4,'LB'],[31.4,31.5,'LC'],[31.5,31.6,'LD'],[31.6,31.7,'LE'],[31.7,31.8,'LF'],[31.8,31.9,'LG'],[31.9,32,'LH'],[32,32.1,'LI'],[32.1,32.2,'LJ'],[32.2,32.3,'LK'],[32.3,32.4,'LL'],[32.4,32.5,'LM'],[32.5,32.6,'LN'],[32.6,32.7,'LO'],[32.7,32.8,'LP'],[32.8,32.9,'LQ'],[32.9,33,'LR'],[33,33.1,'LS'],[33.1,33.2,'LT'],[33.2,33.3,'LU'],[33.3,33.4,'LV'],[33.4,33.5,'LW'],[33.5,33.6,'LX'],[33.6,33.7,'LY'],[33.7,33.8,'LZ'],[33.8,33.9,'MA'],[33.9,34,'MB'],[34,34.1,'MC'],[34.1,34.2,'MD'],[34.2,34.3,'ME'],[34.3,34.4,'MF'],[34.4,34.5,'MG'],[34.5,34.6,'MH'],[34.6,34.7,'MI'],[34.7,34.8,'MJ'],[34.8,34.9,'MK'],[34.9,35,'ML'],[35,35.1,'MM'],[35.1,35.2,'MN'],[35.2,35.3,'MO'],[35.3,35.4,'MP'],[35.4,35.5,'MQ'],[35.5,35.6,'MR'],[35.6,35.7,'MS'],[35.7,35.8,'MT'],[35.8,35.9,'MU'],[35.9,36,'MV'],[36,36.1,'MW'],[36.1,36.2,'MX'],[36.2,36.3,'MY'],[36.3,36.4,'MZ'],[36.4,36.5,'NA'],[36.5,36.6,'NB'],[36.6,36.7,'NC'],[36.7,36.8,'ND'],[36.8,36.9,'NE'],[36.9,37,'NF'],[37,37.1,'NG'],[37.1,37.2,'NH'],[37.2,37.3,'NI'],[37.3,37.4,'NJ'],[37.4,37.5,'NK'],[37.5,37.6,'NL'],[37.6,37.7,'NM'],[37.7,37.8,'NN'],[37.8,37.9,'NO'],[37.9,38,'NP'],[38,38.1,'NQ'],[38.1,38.2,'NR'],[38.2,38.3,'NS'],[38.3,38.4,'NT'],[38.4,38.5,'NU'],[38.5,38.6,'NV'],[38.6,38.7,'NW'],[38.7,38.8,'NX'],[38.8,38.9,'NY'],[38.9,39,'NZ'],[39,39.1,'OA'],[39.1,39.2,'OB'],[39.2,39.3,'OC'],[39.3,39.4,'OD'],[39.4,39.5,'OE'],[39.5,39.6,'OF'],[39.6,39.7,'OG'],[39.7,39.8,'OH'],[39.8,39.9,'OI'],[39.9,40,'OJ'],[40,40.1,'OK'],[40.1,40.2,'OL'],[40.2,40.3,'OM'],[40.3,40.4,'ON'],[40.4,40.5,'OO'],[40.5,40.6,'OP'],[40.6,40.7,'OQ'],[40.7,40.8,'OR'],[40.8,40.9,'OS'],[40.9,41,'OT'],[41,41.1,'OU'],[41.1,41.2,'OV'],[41.2,41.3,'OW'],[41.3,41.4,'OX'],[41.4,41.5,'OY'],[41.5,41.6,'OZ'],[41.6,41.7,'PA'],[41.7,41.8,'PB'],[41.8,41.9,'PC'],[41.9,42,'PD'],[42,42.1,'PE'],[42.1,42.2,'PF'],[42.2,42.3,'PG'],[42.3,42.4,'PH'],[42.4,42.5,'PI'],[42.5,42.6,'PJ'],[42.6,42.7,'PK'],[42.7,42.8,'PL'],[42.8,42.9,'PM'],[42.9,43,'PN'],[43,43.1,'PO'],[43.1,43.2,'PP'],[43.2,43.3,'PQ'],[43.3,43.4,'PR'],[43.4,43.5,'PS'],[43.5,43.6,'PT'],[43.6,43.7,'PU'],[43.7,43.8,'PV'],[43.8,43.9,'PW'],[43.9,44,'PX'],[44,44.1,'PY'],[44.1,44.2,'PZ'],[44.2,44.3,'QA'],[44.3,44.4,'QB'],[44.4,44.5,'QC'],[44.5,44.6,'QD'],[44.6,44.7,'QE'],[44.7,44.8,'QF'],[44.8,44.9,'QG'],[44.9,45,'QH'],[45,45.1,'QI'],[45.1,45.2,'QJ'],[45.2,45.3,'QK'],[45.3,45.4,'QL'],[45.4,45.5,'QM'],[45.5,45.6,'QN'],[45.6,45.7,'QO'],[45.7,45.8,'QP'],[45.8,45.9,'QQ'],[45.9,46,'QR'],[46,46.1,'QS'],[46.1,46.2,'QT'],[46.2,46.3,'QU'],[46.3,46.4,'QV'],[46.4,46.5,'QW'],[46.5,46.6,'QX'],[46.6,46.7,'QY'],[46.7,46.8,'QZ'],[46.8,46.9,'RA'],[46.9,47,'RB'],[47,47.1,'RC'],[47.1,47.2,'RD'],[47.2,47.3,'RE'],[47.3,47.4,'RF'],[47.4,47.5,'RG'],[47.5,47.6,'RH'],[47.6,47.7,'RI'],[47.7,47.8,'RJ'],[47.8,47.9,'RK'],[47.9,48,'RL'],[48,48.1,'RM'],[48.1,48.2,'RN'],[48.2,48.3,'RO'],[48.3,48.4,'RP'],[48.4,48.5,'RQ'],[48.5,48.6,'RR'],[48.6,48.7,'RS'],[48.7,48.8,'RT'],[48.8,48.9,'RU'],[48.9,49,'RV'],[49,49.1,'RW'],[49.1,49.2,'RX'],[49.2,49.3,'RY'],[49.3,49.4,'RZ'],[49.4,49.5,'SA'],[49.5,49.6,'SB'],[49.6,49.7,'SC'],[49.7,49.8,'SD'],[49.8,49.9,'SE'],[49.9,50,'SF'],[50,-0.1,'SG'],[-0.1,-0.2,'SH'],[-0.2,-0.3,'SI'],[-0.3,-0.4,'SJ'],[-0.4,-0.5,'SK'],[-0.5,-0.6,'SL'],[-0.6,-0.7,'SM'],[-0.7,-0.8,'SN'],[-0.8,-0.9,'SO'],[-0.9,-1,'SP'],[-1,-1.1,'SQ'],[-1.1,-1.2,'SR'],[-1.2,-1.3,'SS'],[-1.3,-1.4,'ST'],[-1.4,-1.5,'SU'],[-1.5,-1.6,'SV'],[-1.6,-1.7,'SW'],[-1.7,-1.8,'SX'],[-1.8,-1.9,'SY'],[-1.9,-2,'SZ'],[-2,-2.1,'TA'],[-2.1,-2.2,'TB'],[-2.2,-2.3,'TC'],[-2.3,-2.4,'TD'],[-2.4,-2.5,'TE'],[-2.5,-2.6,'TF'],[-2.6,-2.7,'TG'],[-2.7,-2.8,'TH'],[-2.8,-2.9,'TI'],[-2.9,-3,'TJ'],[-3,-3.1,'TK'],[-3.1,-3.2,'TL'],[-3.2,-3.3,'TM'],[-3.3,-3.4,'TN'],[-3.4,-3.5,'TO'],[-3.5,-3.6,'TP'],[-3.6,-3.7,'TQ'],[-3.7,-3.8,'TR'],[-3.8,-3.9,'TS'],[-3.9,-4,'TT'],[-4,-4.1,'TU'],[-4.1,-4.2,'TV'],[-4.2,-4.3,'TW'],[-4.3,-4.4,'TX'],[-4.4,-4.5,'TY'],[-4.5,-4.6,'TZ'],[-4.6,-4.7,'UA'],[-4.7,-4.8,'UB'],[-4.8,-4.9,'UC'],[-4.9,-5,'UD'],[-5,-5.1,'UE'],[-5.1,-5.2,'UF'],[-5.2,-5.3,'UG'],[-5.3,-5.4,'UH'],[-5.4,-5.5,'UI'],[-5.5,-5.6,'UJ'],[-5.6,-5.7,'UK'],[-5.7,-5.8,'UL'],[-5.8,-5.9,'UM'],[-5.9,-6,'UN'],[-6,-6.1,'UO'],[-6.1,-6.2,'UP'],[-6.2,-6.3,'UQ'],[-6.3,-6.4,'UR'],[-6.4,-6.5,'US'],[-6.5,-6.6,'UT'],[-6.6,-6.7,'UU'],[-6.7,-6.8,'UV'],[-6.8,-6.9,'UW'],[-6.9,-7,'UX'],[-7,-7.1,'UY'],[-7.1,-7.2,'UZ'],[-7.2,-7.3,'VA'],[-7.3,-7.4,'VB'],[-7.4,-7.5,'VC'],[-7.5,-7.6,'VD'],[-7.6,-7.7,'VE'],[-7.7,-7.8,'VF'],[-7.8,-7.9,'VG'],[-7.9,-8,'VH'],[-8,-8.1,'VI'],[-8.1,-8.2,'VJ'],[-8.2,-8.3,'VK'],[-8.3,-8.4,'VL'],[-8.4,-8.5,'VM'],[-8.5,-8.6,'VN'],[-8.6,-8.7,'VO'],[-8.7,-8.8,'VP'],[-8.8,-8.9,'VQ'],[-8.9,-9,'VR'],[-9,-9.1,'VS'],[-9.1,-9.2,'VT'],[-9.2,-9.3,'VU'],[-9.3,-9.4,'VV'],[-9.4,-9.5,'VW'],[-9.5,-9.6,'VX'],[-9.6,-9.7,'VY'],[-9.7,-9.8,'VZ'],[-9.8,-9.9,'WA'],[-9.9,-10,'WB'],[-10,-10.1,'WC'],[-10.1,-10.2,'WD'],[-10.2,-10.3,'WE'],[-10.3,-10.4,'WF'],[-10.4,-10.5,'WG'],[-10.5,-10.6,'WH'],[-10.6,-10.7,'WI'],[-10.7,-10.8,'WJ'],[-10.8,-10.9,'WK'],[-10.9,-11,'WL'],[-11,-11.1,'WM'],[-11.1,-11.2,'WN'],[-11.2,-11.3,'WO'],[-11.3,-11.4,'WP'],[-11.4,-11.5,'WQ'],[-11.5,-11.6,'WR'],[-11.6,-11.7,'WS'],[-11.7,-11.8,'WT'],[-11.8,-11.9,'WU'],[-11.9,-12,'WV'],[-12,-12.1,'WW'],[-12.1,-12.2,'WX'],[-12.2,-12.3,'WY'],[-12.3,-12.4,'WZ'],[-12.4,-12.5,'XA'],[-12.5,-12.6,'XB'],[-12.6,-12.7,'XC'],[-12.7,-12.8,'XD'],[-12.8,-12.9,'XE'],[-12.9,-13,'XF'],[-13,-13.1,'XG'],[-13.1,-13.2,'XH'],[-13.2,-13.3,'XI'],[-13.3,-13.4,'XJ'],[-13.4,-13.5,'XK'],[-13.5,-13.6,'XL'],[-13.6,-13.7,'XM'],[-13.7,-13.8,'XN'],[-13.8,-13.9,'XO'],[-13.9,-14,'XP'],[-14,-14.1,'XQ'],[-14.1,-14.2,'XR'],[-14.2,-14.3,'XS'],[-14.3,-14.4,'XT'],[-14.4,-14.5,'XU'],[-14.5,-14.6,'XV'],[-14.6,-14.7,'XW'],[-14.7,-14.8,'XX'],[-14.8,-14.9,'XY'],[-14.9,-15,'XZ'],[-15,-15.1,'YA'],[-15.1,-15.2,'YB'],[-15.2,-15.3,'YC'],[-15.3,-15.4,'YD'],[-15.4,-15.5,'YE'],[-15.5,-15.6,'YF'],[-15.6,-15.7,'YG'],[-15.7,-15.8,'YH'],[-15.8,-15.9,'YI'],[-15.9,-16,'YJ'],[-16,-16.1,'YK'],[-16.1,-16.2,'YL'],[-16.2,-16.3,'YM'],[-16.3,-16.4,'YN'],[-16.4,-16.5,'YO'],[-16.5,-16.6,'YP'],[-16.6,-16.7,'YQ'],[-16.7,-16.8,'YR'],[-16.8,-16.9,'YS'],[-16.9,-17,'YT'],[-17,-17.1,'YU'],[-17.1,-17.2,'YV'],[-17.2,-17.3,'YW'],[-17.3,-17.4,'YX'],[-17.4,-17.5,'YY'],[-17.5,-17.6,'YZ'],[-17.6,-17.7,'ZA'],[-17.7,-17.8,'ZB'],[-17.8,-17.9,'ZC'],[-17.9,-18,'ZD'],[-18,-18.1,'ZE'],[-18.1,-18.2,'ZF'],[-18.2,-18.3,'ZG'],[-18.3,-18.4,'ZH'],[-18.4,-18.5,'ZI'],[-18.5,-18.6,'ZJ'],[-18.6,-18.7,'ZK'],[-18.7,-18.8,'ZL'],[-18.8,-18.9,'ZM'],[-18.9,-19,'ZN'],[-19,-19.1,'ZO'],[-19.1,-19.2,'ZP'],[-19.2,-19.3,'ZQ'],[-19.3,-19.4,'ZR'],[-19.4,-19.5,'ZS'],[-19.5,-19.6,'ZT'],[-19.6,-19.7,'ZU'],[-19.7,-19.8,'ZV'],[-19.8,-19.9,'ZW'],[-19.9,-20,'ZX'],[-20,-20.1,'ZY'],[-20.1,-20.2,'ZZ'],[-20.2,-20.3,'AAA'],[-20.3,-20.4,'AAB'],[-20.4,-20.5,'AAC'],[-20.5,-20.6,'AAD'],[-20.6,-20.7,'AAE'],[-20.7,-20.8,'AAF'],[-20.8,-20.9,'AAG'],[-20.9,-21,'AAH'],[-21,-21.1,'AAI'],[-21.1,-21.2,'AAJ'],[-21.2,-21.3,'AAK'],[-21.3,-21.4,'AAL'],[-21.4,-21.5,'AAM'],[-21.5,-21.6,'AAN'],[-21.6,-21.7,'AAO'],[-21.7,-21.8,'AAP'],[-21.8,-21.9,'AAQ'],[-21.9,-22,'AAR'],[-22,-22.1,'AAS'],[-22.1,-22.2,'AAT'],[-22.2,-22.3,'AAU'],[-22.3,-22.4,'AAV'],[-22.4,-22.5,'AAW'],[-22.5,-22.6,'AAX'],[-22.6,-22.7,'AAY'],[-22.7,-22.8,'AAZ'],[-22.8,-22.9,'ABA'],[-22.9,-23,'ABB'],[-23,-23.1,'ABC'],[-23.1,-23.2,'ABD'],[-23.2,-23.3,'ABE'],[-23.3,-23.4,'ABF'],[-23.4,-23.5,'ABG'],[-23.5,-23.6,'ABH'],[-23.6,-23.7,'ABI'],[-23.7,-23.8,'ABJ'],[-23.8,-23.9,'ABK'],[-23.9,-24,'ABL'],[-24,-24.1,'ABM'],[-24.1,-24.2,'ABN'],[-24.2,-24.3,'ABO'],[-24.3,-24.4,'ABP'],[-24.4,-24.5,'ABQ'],[-24.5,-24.6,'ABR'],[-24.6,-24.7,'ABS'],[-24.7,-24.8,'ABT'],[-24.8,-24.9,'ABU'],[-24.9,-25,'ABV'],[-25,-25.1,'ABW'],[-25.1,-25.2,'ABX'],[-25.2,-25.3,'ABY'],[-25.3,-25.4,'ABZ'],[-25.4,-25.5,'ACA'],[-25.5,-25.6,'ACB'],[-25.6,-25.7,'ACC'],[-25.7,-25.8,'ACD'],[-25.8,-25.9,'ACE'],[-25.9,-26,'ACF'],[-26,-26.1,'ACG'],[-26.1,-26.2,'ACH'],[-26.2,-26.3,'ACI'],[-26.3,-26.4,'ACJ'],[-26.4,-26.5,'ACK'],[-26.5,-26.6,'ACL'],[-26.6,-26.7,'ACM'],[-26.7,-26.8,'ACN'],[-26.8,-26.9,'ACO'],[-26.9,-27,'ACP'],[-27,-27.1,'ACQ'],[-27.1,-27.2,'ACR'],[-27.2,-27.3,'ACS'],[-27.3,-27.4,'ACT'],[-27.4,-27.5,'ACU'],[-27.5,-27.6,'ACV'],[-27.6,-27.7,'ACW'],[-27.7,-27.8,'ACX'],[-27.8,-27.9,'ACY'],[-27.9,-28,'ACZ'],[-28,-28.1,'ADA'],[-28.1,-28.2,'ADB'],[-28.2,-28.3,'ADC'],[-28.3,-28.4,'ADD'],[-28.4,-28.5,'ADE'],[-28.5,-28.6,'ADF'],[-28.6,-28.7,'ADG'],[-28.7,-28.8,'ADH'],[-28.8,-28.9,'ADI'],[-28.9,-29,'ADJ'],[-29,-29.1,'ADK'],[-29.1,-29.2,'ADL'],[-29.2,-29.3,'ADM'],[-29.3,-29.4,'ADN'],[-29.4,-29.5,'ADO'],[-29.5,-29.6,'ADP'],[-29.6,-29.7,'ADQ'],[-29.7,-29.8,'ADR'],[-29.8,-29.9,'ADS'],[-29.9,-30,'ADT'],[-30,-30.1,'ADU'],[-30.1,-30.2,'ADV'],[-30.2,-30.3,'ADW'],[-30.3,-30.4,'ADX'],[-30.4,-30.5,'ADY'],[-30.5,-30.6,'ADZ'],[-30.6,-30.7,'AEA'],[-30.7,-30.8,'AEB'],[-30.8,-30.9,'AEC'],[-30.9,-31,'AED'],[-31,-31.1,'AEE'],[-31.1,-31.2,'AEF'],[-31.2,-31.3,'AEG'],[-31.3,-31.3999999999999,'AEH'],[-31.3999999999999,-31.4999999999999,'AEI'],[-31.4999999999999,-31.5999999999999,'AEJ'],[-31.5999999999999,-31.6999999999999,'AEK'],[-31.6999999999999,-31.7999999999999,'AEL'],[-31.7999999999999,-31.8999999999999,'AEM'],[-31.8999999999999,-31.9999999999999,'AEN'],[-31.9999999999999,-32.0999999999999,'AEO'],[-32.0999999999999,-32.1999999999999,'AEP'],[-32.1999999999999,-32.2999999999999,'AEQ'],[-32.2999999999999,-32.3999999999999,'AER'],[-32.3999999999999,-32.4999999999999,'AES'],[-32.4999999999999,-32.5999999999999,'AET'],[-32.5999999999999,-32.6999999999999,'AEU'],[-32.6999999999999,-32.7999999999999,'AEV'],[-32.7999999999999,-32.8999999999999,'AEW'],[-32.8999999999999,-32.9999999999999,'AEX'],[-32.9999999999999,-33.0999999999999,'AEY'],[-33.0999999999999,-33.1999999999999,'AEZ'],[-33.1999999999999,-33.2999999999999,'AFA'],[-33.2999999999999,-33.3999999999999,'AFB'],[-33.3999999999999,-33.4999999999999,'AFC'],[-33.4999999999999,-33.5999999999999,'AFD'],[-33.5999999999999,-33.6999999999999,'AFE'],[-33.6999999999999,-33.7999999999999,'AFF'],[-33.7999999999999,-33.8999999999999,'AFG'],[-33.8999999999999,-33.9999999999999,'AFH'],[-33.9999999999999,-34.0999999999999,'AFI'],[-34.0999999999999,-34.1999999999999,'AFJ'],[-34.1999999999999,-34.2999999999999,'AFK'],[-34.2999999999999,-34.3999999999999,'AFL'],[-34.3999999999999,-34.4999999999999,'AFM'],[-34.4999999999999,-34.5999999999999,'AFN'],[-34.5999999999999,-34.6999999999999,'AFO'],[-34.6999999999999,-34.7999999999999,'AFP'],[-34.7999999999999,-34.8999999999999,'AFQ'],[-34.8999999999999,-34.9999999999999,'AFR'],[-34.9999999999999,-35.0999999999999,'AFS'],[-35.0999999999999,-35.1999999999999,'AFT'],[-35.1999999999999,-35.2999999999999,'AFU'],[-35.2999999999999,-35.3999999999999,'AFV'],[-35.3999999999999,-35.4999999999999,'AFW'],[-35.4999999999999,-35.5999999999999,'AFX'],[-35.5999999999999,-35.6999999999999,'AFY'],[-35.6999999999999,-35.7999999999999,'AFZ'],[-35.7999999999999,-35.8999999999999,'AGA'],[-35.8999999999999,-35.9999999999999,'AGB'],[-35.9999999999999,-36.0999999999999,'AGC'],[-36.0999999999999,-36.1999999999999,'AGD'],[-36.1999999999999,-36.2999999999999,'AGE'],[-36.2999999999999,-36.3999999999999,'AGF'],[-36.3999999999999,-36.4999999999999,'AGG'],[-36.4999999999999,-36.5999999999999,'AGH'],[-36.5999999999999,-36.6999999999999,'AGI'],[-36.6999999999999,-36.7999999999999,'AGJ'],[-36.7999999999999,-36.8999999999999,'AGK'],[-36.8999999999999,-36.9999999999999,'AGL'],[-36.9999999999999,-37.0999999999999,'AGM'],[-37.0999999999999,-37.1999999999999,'AGN'],[-37.1999999999999,-37.2999999999999,'AGO'],[-37.2999999999999,-37.3999999999999,'AGP'],[-37.3999999999999,-37.4999999999999,'AGQ'],[-37.4999999999999,-37.5999999999999,'AGR'],[-37.5999999999999,-37.6999999999999,'AGS'],[-37.6999999999999,-37.7999999999999,'AGT'],[-37.7999999999999,-37.8999999999999,'AGU'],[-37.8999999999999,-37.9999999999999,'AGV'],[-37.9999999999999,-38.0999999999999,'AGW'],[-38.0999999999999,-38.1999999999999,'AGX'],[-38.1999999999999,-38.2999999999999,'AGY'],[-38.2999999999999,-38.3999999999999,'AGZ'],[-38.3999999999999,-38.4999999999999,'AHA'],[-38.4999999999999,-38.5999999999999,'AHB'],[-38.5999999999999,-38.6999999999999,'AHC'],[-38.6999999999999,-38.7999999999999,'AHD'],[-38.7999999999999,-38.8999999999999,'AHE'],[-38.8999999999999,-38.9999999999999,'AHF'],[-38.9999999999999,-39.0999999999999,'AHG'],[-39.0999999999999,-39.1999999999999,'AHH'],[-39.1999999999999,-39.2999999999999,'AHI'],[-39.2999999999999,-39.3999999999999,'AHJ'],[-39.3999999999999,-39.4999999999999,'AHK'],[-39.4999999999999,-39.5999999999999,'AHL'],[-39.5999999999999,-39.6999999999999,'AHM'],[-39.6999999999999,-39.7999999999999,'AHN'],[-39.7999999999999,-39.8999999999999,'AHO'],[-39.8999999999999,-39.9999999999999,'AHP'],[-39.9999999999999,-40.0999999999999,'AHQ'],[-40.0999999999999,-40.1999999999999,'AHR'],[-40.1999999999999,-40.2999999999999,'AHS'],[-40.2999999999999,-40.3999999999999,'AHT'],[-40.3999999999999,-40.4999999999999,'AHU'],[-40.4999999999999,-40.5999999999999,'AHV'],[-40.5999999999999,-40.6999999999999,'AHW'],[-40.6999999999999,-40.7999999999999,'AHX'],[-40.7999999999999,-40.8999999999999,'AHY'],[-40.8999999999999,-40.9999999999999,'AHZ'],[-40.9999999999999,-41.0999999999999,'AIA'],[-41.0999999999999,-41.1999999999999,'AIB'],[-41.1999999999999,-41.2999999999999,'AIC'],[-41.2999999999999,-41.3999999999999,'AID'],[-41.3999999999999,-41.4999999999999,'AIE'],[-41.4999999999999,-41.5999999999999,'AIF'],[-41.5999999999999,-41.6999999999999,'AIG'],[-41.6999999999999,-41.7999999999999,'AIH'],[-41.7999999999999,-41.9,'AII'],[-41.9,-42,'AIJ'],[-42,-42.1,'AIK'],[-42.1,-42.2,'AIL'],[-42.2,-42.3,'AIM'],[-42.3,-42.4,'AIN'],[-42.4,-42.5,'AIO'],[-42.5,-42.6,'AIP'],[-42.6,-42.7,'AIQ'],[-42.7,-42.8,'AIR'],[-42.8,-42.9,'AIS'],[-42.9,-43,'AIT'],[-43,-43.1,'AIU'],[-43.1,-43.2,'AIV'],[-43.2,-43.3,'AIW'],[-43.3,-43.4,'AIX'],[-43.4,-43.5,'AIY'],[-43.5,-43.6,'AIZ'],[-43.6,-43.7,'AJA'],[-43.7,-43.8,'AJB'],[-43.8,-43.9,'AJC'],[-43.9,-44,'AJD'],[-44,-44.1,'AJE'],[-44.1,-44.2,'AJF'],[-44.2,-44.3,'AJG'],[-44.3,-44.4,'AJH'],[-44.4,-44.5,'AJI'],[-44.5,-44.6,'AJJ'],[-44.6,-44.7,'AJK'],[-44.7,-44.8,'AJL'],[-44.8,-44.9,'AJM'],[-44.9,-45,'AJN'],[-45,-45.1,'AJO'],[-45.1,-45.2,'AJP'],[-45.2,-45.3,'AJQ'],[-45.3,-45.4,'AJR'],[-45.4,-45.5,'AJS'],[-45.5,-45.6,'AJT'],[-45.6,-45.7,'AJU'],[-45.7,-45.8,'AJV'],[-45.8,-45.9,'AJW'],[-45.9,-46,'AJX'],[-46,-46.1,'AJY'],[-46.1,-46.2,'AJZ'],[-46.2,-46.3,'AKA'],[-46.3,-46.4,'AKB'],[-46.4,-46.5,'AKC'],[-46.5,-46.6,'AKD'],[-46.6,-46.7,'AKE'],[-46.7,-46.8,'AKF'],[-46.8,-46.9,'AKG'],[-46.9,-47,'AKH'],[-47,-47.1,'AKI'],[-47.1,-47.2,'AKJ'],[-47.2,-47.3,'AKK'],[-47.3,-47.4,'AKL'],[-47.4,-47.5,'AKM'],[-47.5,-47.6,'AKN'],[-47.6,-47.7,'AKO'],[-47.7,-47.8,'AKP'],[-47.8,-47.9,'AKQ'],[-47.9,-48,'AKR'],[-48,-48.1,'AKS'],[-48.1,-48.2,'AKT'],[-48.2,-48.3,'AKU'],[-48.3,-48.4,'AKV'],[-48.4,-48.5000000000002,'AKW'],[-48.5000000000002,-48.6000000000002,'AKX'],[-48.6000000000002,-48.7000000000002,'AKY'],[-48.7000000000002,-48.8000000000002,'AKZ'],[-48.8000000000002,-48.9000000000002,'ALA'],[-48.9000000000002,-49.0000000000002,'ALB'],[-49.0000000000002,-49.1000000000002,'ALC'],[-49.1000000000002,-49.2000000000002,'ALD'],[-49.2000000000002,-49.3000000000002,'ALE'],[-49.3000000000002,-49.4000000000002,'ALF'],[-49.4000000000002,-49.5000000000002,'ALG'],[-49.5000000000002,-49.6000000000002,'ALH'],[-49.6000000000002,-49.7000000000002,'ALI'],[-49.7000000000002,-49.8000000000002,'ALJ'],[-49.8000000000002,-49.9000000000002,'ALK'],[-49.9000000000002,-50.0000000000002,'ALL']]

path = "/Users/benjaminsamudio/ReLieF_Fingerprints_Demonstration/Abl_1b_6xr6_0006_ATPbindingSite/"
os.chdir(path)

header_string = "Name" + "," + "Color" + "," + "Fingerprints" + ","

# for header_index in range(0,10000): # This assumes that there are 10,000 bits
#    if header_index == 0:
#        header_string = header_string + "Bit_" + str(header_index) + " "
#    elif header_index > 0:
#        header_string = header_string + " " + "Bit_" + str(header_index)

print(header_string)

fingerprints_filename = path + "ReLieF_Fingerprints.csv"
    

with open(fingerprints_filename,'w') as file_object:
    header_string = header_string + "\n"
    file_object.write(header_string)

for file in os.listdir(): # START OF ITERATIONS
    if file.endswith(".wrl"):
        print(f"Currently processing this file: {file}")
        temp_filename_tuple = ()
        temp_filename_tuple = os.path.splitext(file)
        print(f"The filename, without extension, is: {temp_filename_tuple[0]}")
        output_filename_base = temp_filename_tuple[0]
        # test_surface_file = '/Users/benjaminsamudio/ReLieF_Fingerprints_Demonstration/PDB_5r7y_cleaned_reoriented_setView.wrl'
        test_surface_file = file
        test_output_file = path + 'ReLieF_ReceptorSurfaceSegments_' + output_filename_base +  '.xyz'
        output_filename = path + 'ReLieF_ReceptorFingerprintBitLocations_' + output_filename_base + '.xyz'
        structure_fingerprint_correspondence = path + 'ReLieF_StructureFingerprintBitCorrespondence_' + output_filename_base + '.csv'
        map_filename = path + 'ReLieF_Map_' + output_filename_base + '.csv'
        fingerprint_output_string = output_filename_base + "," + "FILL_ME" + ","

        with open(test_output_file,'a') as file_object:
            file_object.write("This is a test of the blanket fingerprints\n")
            file_object.write("\n")
            
        with open(structure_fingerprint_correspondence,'a') as file_object:
            file_object.write("map_row,map_column,bit_index,structure_index\n")
            

        ######################################## Extract the vertex points constituting surface triangles in the *wrl file 

        triangle_vertex_end_flag = 0
        vertex_concatenate_count = 0
        vertex_coordinate_point = []
        vertex_trio = []

        with open(test_surface_file) as file_object:
            for surface_triangle_coordinates in file_object:
                triangle_vertex_coordinates = []
                triangle_vertex_end = []
                triangle_vertex_coordinates = re.search(r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+),",surface_triangle_coordinates)
                triangle_vertex_end  = re.search(r"coordIndex",surface_triangle_coordinates)
                if triangle_vertex_end:
                    triangle_vertex_end_flag = 1
                if triangle_vertex_coordinates and triangle_vertex_end_flag == 0:
                    vertex_concatenate_count += 1
                    vertex_coordinate_point.append(float(triangle_vertex_coordinates.group(1)))
                    vertex_coordinate_point.append(float(triangle_vertex_coordinates.group(2)))
                    vertex_coordinate_point.append(float(triangle_vertex_coordinates.group(3)))
                #   print(vertex_line_string.group(1),vertex_line_string.group(2),vertex_line_string.group(3))
                    if vertex_concatenate_count == 3:
                        vertex_trio.append(vertex_coordinate_point)
                #       print(test_count,vertex_points_string)
                        vertex_concatenate_count = 0;
                        vertex_coordinate_point = []

        ######################################## Determine the intersection points, if any, between segment planes and triangle vectors 

        main_normal_x_start = 0
        main_normal_y_start = 0
        main_normal_z_start = 0

        main_normal_x_end = 1
        main_normal_y_end = 0
        main_normal_z_end = 0

        main_normal_x = main_normal_x_end - main_normal_x_start # This is variable "A".  The normal vector specifies the direction.  In this case, along the x-axis
        main_normal_y = main_normal_y_end - main_normal_y_start # This is variable "B".  The normal vector specifies the direction.  In this case, along the x-axis
        main_normal_z = main_normal_z_end - main_normal_z_start # This is variable "C".  The normal vector specifies the direction.  In this case, along the x-axis

        main_x_start = -50
        main_y_start = 0
        main_z_start = 0

        main_x_end = 50
        main_y_end = 0
        main_z_end = 0

        # triangle_vector_x_start = 0
        # triangle_vector_y_start = 0
        # triangle_vector_z_start = 0

        # triangle_vector_x_end = 10
        # triangle_vector_y_end = 0
        # triangle_vector_z_end = 0 

        # triangle_vector_x = triangle_vector_x_end - triangle_vector_x_start 
        # triangle_vector_y = triangle_vector_y_end - triangle_vector_y_start
        # triangle_vector_z = triangle_vector_z_end - triangle_vector_z_start

        total_blanket_rows = 100 # This is the total number of rows in the "blanket".  Coordinate this value with the main x, y, and z coordinates
        intersect_points = []

        for blanket_row_index in range(0,total_blanket_rows): # Create the "blanket" rows
            blanket_row = main_normal_x * (main_x_start + (main_x_end - main_x_start) * blanket_row_index/total_blanket_rows) \
            + main_normal_y * (main_y_start + (main_y_end - main_y_start) * blanket_row_index/total_blanket_rows) 
            + main_normal_z * (main_z_start + (main_z_end - main_z_start) * blanket_row_index/total_blanket_rows) # This is variable "D".
            intersect_points_components = []
        #   intersect_points = [] # Disabled on 04FEB2023 so that all intersection points could be collected into one list
            for triangle_vertex_trio in vertex_trio:
                for vertex_trio_index in range(0,3):
                    if vertex_trio_index == 0:
                        triangle_vector_x_start = triangle_vertex_trio[0]
                        triangle_vector_y_start = triangle_vertex_trio[1]
                        triangle_vector_z_start = triangle_vertex_trio[2]
                        triangle_vector_x = triangle_vertex_trio[3] - triangle_vertex_trio[0]  
                        triangle_vector_y = triangle_vertex_trio[4] - triangle_vertex_trio[1]
                        triangle_vector_z = triangle_vertex_trio[5] - triangle_vertex_trio[2]
                    elif vertex_trio_index == 1:
                        triangle_vector_x_start = triangle_vertex_trio[0]
                        triangle_vector_y_start = triangle_vertex_trio[1]
                        triangle_vector_z_start = triangle_vertex_trio[2]
                        triangle_vector_x = triangle_vertex_trio[6] - triangle_vertex_trio[0]  
                        triangle_vector_y = triangle_vertex_trio[7] - triangle_vertex_trio[1]
                        triangle_vector_z = triangle_vertex_trio[8] - triangle_vertex_trio[2]
                    elif vertex_trio_index == 2:
                        triangle_vector_x_start = triangle_vertex_trio[3]
                        triangle_vector_y_start = triangle_vertex_trio[4]
                        triangle_vector_z_start = triangle_vertex_trio[5]
                        triangle_vector_x = triangle_vertex_trio[6] - triangle_vertex_trio[3]  
                        triangle_vector_y = triangle_vertex_trio[7] - triangle_vertex_trio[4]
                        triangle_vector_z = triangle_vertex_trio[8] - triangle_vertex_trio[5]
                    if (main_normal_x * triangle_vector_x  + main_normal_y * triangle_vector_y  + main_normal_z * triangle_vector_z) != 0:
                        t_parameter = -1 * (main_normal_x * triangle_vector_x_start + main_normal_y * triangle_vector_y_start + main_normal_z * triangle_vector_z_start + blanket_row) \
                        / (main_normal_x * triangle_vector_x  + main_normal_y * triangle_vector_y  + main_normal_z * triangle_vector_z)
                    #   print(blanket_row,t_parameter)
                        if abs(t_parameter) < 1 and abs(t_parameter) > 0 or t_parameter == 1 or t_parameter == 0 or t_parameter == -1:
                            intersect_point_x = triangle_vector_x_start - (triangle_vector_x * (main_normal_x * triangle_vector_x_start + main_normal_y * triangle_vector_y_start + main_normal_z * triangle_vector_z_start + blanket_row)/(main_normal_x * triangle_vector_x  + main_normal_y * triangle_vector_y  + main_normal_z * triangle_vector_z))
                            intersect_point_y = triangle_vector_y_start - (triangle_vector_y * (main_normal_x * triangle_vector_x_start + main_normal_y * triangle_vector_y_start + main_normal_z * triangle_vector_z_start + blanket_row)/(main_normal_x * triangle_vector_x  + main_normal_y * triangle_vector_y  + main_normal_z * triangle_vector_z))
                            intersect_point_z = triangle_vector_z_start - (triangle_vector_z * (main_normal_x * triangle_vector_x_start + main_normal_y * triangle_vector_y_start + main_normal_z * triangle_vector_z_start + blanket_row)/(main_normal_x * triangle_vector_x  + main_normal_y * triangle_vector_y  + main_normal_z * triangle_vector_z))
                            intersect_points_components = [intersect_point_x,intersect_point_y,intersect_point_z]
                            intersect_points.append(intersect_points_components)
                    #       print(intersect_points[0])
                    #       print(blanket_row,t_parameter,intersect_point_x,intersect_point_y,intersect_point_z)
                            with open(test_output_file,'a') as file_object:
                                file_object.write(f"Cl {intersect_point_x} {intersect_point_y} {intersect_point_z}\n")
                        else:
                            next



        y_axis_start = -50
        y_axis_end = 50
        z_axis_start = -50
        z_axis_end = 50
        z_axis_width = (abs(z_axis_start)+abs(z_axis_end))
        z_axis_increment_size = 0.5 
        z_axis_max_increments = round(z_axis_width/z_axis_increment_size)
        y_axis_increment = 1
        y_axis_half_increment = y_axis_increment / 2
        sampling_point_coordinates = []
        z_coordinates = []
        sampling_circle_radius = 1
        in_circle_flag = 0
        distance_wiggle = 0.50
        temporary_append = []
        output_row_concatenate = []
        output_bit_concatenate = []
        output_row_count = 0 # X-axis value
        output_column_count = 0 # Y-axis value
        output_signal_count = 0 # A count of one is made if a bit has a value in it (not equal to "#" which is considered NULL)
        fingerprint_bit_count = 0
        
        for x_axis_coordinate in range(main_x_start,main_x_end+1): # Loop through segments
            output_row_concatenate = ['#'] * (abs(y_axis_start) + abs(y_axis_end))
            print(f"Now processing the following segment: {x_axis_coordinate}")
            bit_count = 0
            output_row_count += 1
            output_column_count = 0
            for y_axis_coordinate in range(y_axis_start,y_axis_end): # Move from "right" to "left" of the segments
                old_distance_average = 0
                new_distance_average = 0
                z_coordinates = []
                output_bit_concatenate = []
                intersect_points_y_window = []
                bit_count += 1
                output_column_count += 1
                fingerprint_bit_count += 1
        ###### Extract only those points that are at x_axis_coordinate and within y-axis bounds
                for cluster_candidate_index in range(0,len(intersect_points)): # Loop through all intersection points
                    if intersect_points[cluster_candidate_index][0] == x_axis_coordinate and intersect_points[cluster_candidate_index][1] > (y_axis_coordinate - y_axis_half_increment) and intersect_points[cluster_candidate_index][1] < (y_axis_coordinate + y_axis_half_increment): # Check that the intersection points are available and within the y-axis bounds
                        temporary_append = []
                        temporary_append = [intersect_points[cluster_candidate_index][0],intersect_points[cluster_candidate_index][1],intersect_points[cluster_candidate_index][2]]
                        intersect_points_y_window.append(temporary_append)
        ###### Move a circle from "top" to "bottom" within the y-axis bounds and find the average distance of points in that circle
                for z_axis_increment_index in range(0,z_axis_max_increments): # Move from "top" to "bottom" of the segments
                    z_axis_coordinate = z_axis_start + z_axis_increment_size * z_axis_increment_index
                    distance_average = 0
                    average_distance_count = 0
                    sampling_intersection_distance = 0
                    in_circle_flag = 0
                    for y_window_index in range(0,len(intersect_points_y_window)):
                        sampling_point_coordinates = []
                        intersection_point_coordinates = []
                        sampling_point_coordinates = [x_axis_coordinate,y_axis_coordinate,z_axis_coordinate]
                        intersection_point_coordinates = [intersect_points_y_window[y_window_index][0],intersect_points_y_window[y_window_index][1],intersect_points_y_window[y_window_index][2]]
                        sampling_intersection_distance = math.dist(sampling_point_coordinates,intersection_point_coordinates)
                        if sampling_intersection_distance <= sampling_circle_radius:
                            distance_average = ((distance_average * average_distance_count) + sampling_intersection_distance) / (average_distance_count + 1)
                            average_distance_count += 1
                            in_circle_flag = 1
        ###### Keep only those points constituting an average distance which is within a narrow band within the y-axis bounds
                    if in_circle_flag == 1 and distance_average < distance_wiggle:
                        temporary_append = []
                        temporary_append = [z_axis_coordinate,"AVAILABLE"]
                        z_coordinates.append(temporary_append)
        ###### Keep only one point for each time that the circle moves across the segment contour
                for partner_one_index in range(0,len(z_coordinates)):
                    for partner_two_index in range(partner_one_index + 1,len(z_coordinates)):
                        if abs(z_coordinates[partner_one_index][0] - z_coordinates[partner_two_index][0]) <= distance_wiggle and z_coordinates[partner_one_index][1] == "AVAILABLE" and z_coordinates[partner_two_index][1] == "AVAILABLE":
                            average_z = 0
                            average_z = (z_coordinates[partner_one_index][0] + z_coordinates[partner_two_index][0])/2
                            z_coordinates[partner_one_index][1] = "UNAVAILABLE"
                            z_coordinates[partner_two_index][1] = "UNAVAILABLE"
                            for distance_bin_index in range(0,len(distance_codes)):
                                if average_z >= distance_codes[distance_bin_index][0] and average_z < distance_codes[distance_bin_index][1]:
                                    temporary_append_string = ""
                                    temporary_append_string = distance_codes[distance_bin_index][2] + ":"
                                    output_bit_concatenate.append(temporary_append_string)
                                    break
                            output_signal_count += 1
                            with open(output_filename,'a') as file_object:
                                file_object.write(f"Br {x_axis_coordinate} {y_axis_coordinate} {average_z}\n")
                            with open(structure_fingerprint_correspondence,'a') as file_object:
                                file_object.write(f"{output_row_count},{output_column_count},{fingerprint_bit_count},{output_signal_count}\n")
                for no_partner_index in range(0,len(z_coordinates)):
                    if z_coordinates[no_partner_index][1] == "AVAILABLE":
                        for distance_bin_index in range(0,len(distance_codes)):
                            if z_coordinates[no_partner_index][0] >= distance_codes[distance_bin_index][0] and z_coordinates[no_partner_index][0] < distance_codes[distance_bin_index][1]:
                                temporary_append_string = ''
                                temporary_append_string = distance_codes[distance_bin_index][2] + ":"
                                output_bit_concatenate.append(temporary_append_string)
                                break
                        output_signal_count += 1
                        with open(output_filename,'a') as file_object:
                            file_object.write(f"Br {x_axis_coordinate} {y_axis_coordinate} {z_coordinates[no_partner_index][0]}\n")
                        with open(structure_fingerprint_correspondence,'a') as file_object:
                            file_object.write(f"{output_row_count},{output_column_count},{fingerprint_bit_count},{output_signal_count}\n")    
                if output_bit_concatenate:
                    output_row_concatenate[bit_count] = str(''.join(output_bit_concatenate))
            print(",".join(output_row_concatenate))
            fingerprint_output_string = fingerprint_output_string + " " + " ".join(output_row_concatenate)
            with open(map_filename,'a') as file_object:
                file_object.write(",".join(output_row_concatenate))
                file_object.write("\n")
        with open(fingerprints_filename,'a') as file_object:
                file_object.write(fingerprint_output_string)
                file_object.write("\n")


