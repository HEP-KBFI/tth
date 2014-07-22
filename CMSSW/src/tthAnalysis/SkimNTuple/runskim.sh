 #!/bin/sh
./do.exe  FileToSkim/qcd20      >& log_ ; mv log_ log/qcd20     ; mv output.root Skim/qcd20.root
./do.exe  FileToSkim/tth125     >& log_ ; mv log_ log/tth125    ; mv output.root Skim/tth125.root
./do.exe  FileToSkim/dyjets_m50         >& log_ ; mv log_ log/dyjets_m50        ; mv output.root Skim/dyjets_m50.root
./do.exe  FileToSkim/dyjets_m10to50     >& log_ ; mv log_ log/dyjets_m10to50    ; mv output.root Skim/dyjets_m10to50.root
./do.exe  FileToSkim/wjets      >& log_ ; mv log_ log/wjets     ; mv output.root Skim/wjets.root
./do.exe  FileToSkim/ttjets     >& log_ ; mv log_ log/ttjets    ; mv output.root Skim/ttjets.root
./do.exe  FileToSkim/ttwjets    >& log_ ; mv log_ log/ttwjets   ; mv output.root Skim/ttwjets.root
./do.exe  FileToSkim/ttzjets    >& log_ ; mv log_ log/ttzjets   ; mv output.root Skim/ttzjets.root
./do.exe  FileToSkim/t_t        >& log_ ; mv log_ log/t_t       ; mv output.root Skim/t_t.root
./do.exe  FileToSkim/t_tw       >& log_ ; mv log_ log/t_tw      ; mv output.root Skim/t_tw.root
./do.exe  FileToSkim/t_s        >& log_ ; mv log_ log/t_s       ; mv output.root Skim/t_s.root
./do.exe  FileToSkim/tbar_t     >& log_ ; mv log_ log/tbar_t    ; mv output.root Skim/tbar_t.root
./do.exe  FileToSkim/tbar_tw    >& log_ ; mv log_ log/tbar_tw   ; mv output.root Skim/tbar_tw.root
./do.exe  FileToSkim/tbar_s     >& log_ ; mv log_ log/tbar_s    ; mv output.root Skim/tbar_s.root
./do.exe  FileToSkim/ww         >& log_ ; mv log_ log/ww        ; mv output.root Skim/ww.root
./do.exe  FileToSkim/wz         >& log_ ; mv log_ log/wz        ; mv output.root Skim/wz.root
./do.exe  FileToSkim/zz         >& log_ ; mv log_ log/zz        ; mv output.root Skim/zz.root
./do.exe  FileToSkim/singleMuRun2012A    >& log_ ; mv log_ log/singleMuRun2012A       ; mv output.root Skim/singleMuRun2012A.root
./do.exe  FileToSkim/singleMuRun2012B    >& log_ ; mv log_ log/singleMuRun2012B       ; mv output.root Skim/singleMuRun2012B.root
./do.exe  FileToSkim/singleMuRun2012C    >& log_ ; mv log_ log/singleMuRun2012C       ; mv output.root Skim/singleMuRun2012C.root
./do.exe  FileToSkim/singleMuRun2012D    >& log_ ; mv log_ log/singleMuRun2012D       ; mv output.root Skim/singleMuRun2012D.root

