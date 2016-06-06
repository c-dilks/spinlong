// simply prints a list of the triggers

void ListTriggers() {
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * T = new LevelTwo(RD->env);
  for(int x=0; x<T->N; x++) printf("trigger %s\n",(T->Name(x)).Data());
};
