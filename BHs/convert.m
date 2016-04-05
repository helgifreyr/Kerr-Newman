input  = ToString[$CommandLine[[4]]];
dat = ReadList[input];
Export["gnuplot-"<>input,dat,"Table"]
