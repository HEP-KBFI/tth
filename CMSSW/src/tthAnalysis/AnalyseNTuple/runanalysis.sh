 #!/bin/sh
./do.exe  Skim/tth125.root  > log_   ; mv log_ log/tth125    ; mv output.root hist/tth125.root
./do.exe  Skim/dyjets_m10to50.root > log_    ; mv log_ log/dyjets_m10to50    ; mv output.root hist/dyjets_m10to50.root
./do.exe  Skim/dyjets_m50.root     > log_    ; mv log_ log/dyjets_m50        ; mv output.root hist/dyjets_m50.root
./do.exe  Skim/wjets.root   > log_    ; mv log_ log/wjets     ; mv output.root hist/wjets.root
./do.exe  Skim/ttjets.root  > log_    ; mv log_ log/ttjets    ; mv output.root hist/ttjets.root
./do.exe  Skim/ttwjets.root > log_    ; mv log_ log/ttwjets   ; mv output.root hist/ttwjets.root
./do.exe  Skim/ttzjets.root > log_    ; mv log_ log/ttzjets   ; mv output.root hist/ttzjets.root
./do.exe  Skim/t_t.root     > log_    ; mv log_ log/t_t       ; mv output.root hist/t_t.root
./do.exe  Skim/t_tw.root    > log_    ; mv log_ log/t_tw      ; mv output.root hist/t_tw.root
./do.exe  Skim/t_s.root     > log_    ; mv log_ log/t_s       ; mv output.root hist/t_s.root
./do.exe  Skim/tbar_t.root  > log_    ; mv log_ log/tbar_t    ; mv output.root hist/tbar_t.root
./do.exe  Skim/tbar_tw.root > log_    ; mv log_ log/tbar_tw   ; mv output.root hist/tbar_tw.root
./do.exe  Skim/tbar_s.root  > log_    ; mv log_ log/tbar_s    ; mv output.root hist/tbar_s.root
./do.exe  Skim/ww.root      > log_    ; mv log_ log/ww        ; mv output.root hist/ww.root
./do.exe  Skim/wz.root      > log_    ; mv log_ log/wz        ; mv output.root hist/wz.root
./do.exe  Skim/zz.root      > log_    ; mv log_ log/zz        ; mv output.root hist/zz.root
./do.exe  Skim/singleMuRun2012A.root > log_    ; mv log_ log/singleMuRun2012A       ; mv output.root hist/singleMuRun2012A.root
./do.exe  Skim/singleMuRun2012B.root > log_    ; mv log_ log/singleMuRun2012B       ; mv output.root hist/singleMuRun2012B.root
./do.exe  Skim/singleMuRun2012C.root > log_    ; mv log_ log/singleMuRun2012C       ; mv output.root hist/singleMuRun2012C.root
./do.exe  Skim/singleMuRun2012D.root > log_    ; mv log_ log/singleMuRun2012D       ; mv output.root hist/singleMuRun2012D.root
#./do.exe  Skim/qcdA.root > log_    ; mv log_ log/singleMuRun2012Aqcd       ; mv output.root hist/singleMuRun2012Aqcd.root
#./do.exe  Skim/qcdB.root > log_    ; mv log_ log/singleMuRun2012Bqcd       ; mv output.root hist/singleMuRun2012Bqcd.root
#./do.exe  Skim/qcdC.root > log_    ; mv log_ log/singleMuRun2012Cqcd       ; mv output.root hist/singleMuRun2012Cqcd.root
#./do.exe  Skim/qcdD.root > log_    ; mv log_ log/singleMuRun2012Dqcd       ; mv output.root hist/singleMuRun2012Dqcd.root

