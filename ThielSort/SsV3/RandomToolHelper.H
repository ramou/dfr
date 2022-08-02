/* $Header:  $ */ 

#ifndef RandomToolHelperH
#define RandomToolHelperH

/********************************************************************
 *
 *                            RandomToolHelper.H
 *
 * Functions to support creation and configuration of the random
 * generators. Included into namespace RandomTool in RandomTool.H
 ********************************************************************/

// Given a seed spec, return a seed
SeedT GetSeed(unsigned long long t) {
  if (t == 0) return EngineT::default_seed;
  else if (t == 1) {
    std::random_device rd;
    return rd();
  }
  return t;  // Lose high order 0 bits!
}
void PrintStringList(const char**const list_v) {
  std::cout << list_v[0];
  for (int i = 1; list_v[i][0] != 0; i++) std::cout << list_v[i];
}
void PrintDistHelps();
// Print help about cmd line requests for random dist, configuration.
void CmdLineHelper() {
  std::cout <<
    "Generate random numbers using std::mt19937. Command line spec:\n";
  std::cout << "Generator [seed_spec [configuration]] ... \n";
  std::cout << "Where Generator is one of: ";
  std::cout << GenNames[0];
  for (int i = 1; GenNames[i][0] != 0; i++) std::cout << ", " << GenNames[i];
  std::cout << " . Default is " << DefaultGenerator << '\n';
  const char* SeedOptions[] =
    { "seed spec is one of:\n",
      "0 : Use EngineT::default_seed as the seed\n",
      "1  : Use a random seed from std::random_device rd (the default)\n",
      "any other value Will be used as the seed.\n",
      "configuration varies with distribution and contains 0",
      " or more values\n",
      "and ... is the rest of the arg list (0 or more args). \n",
      "If any optional values are omitted, the '...' must be empty or start ",
      "with an argument whose first character is '-'.\n",
      "" };
  PrintStringList(SeedOptions); std::cout << std::endl;
  PrintDistHelps();
}
// return 0 iff *s1 == *s2.
int strcmpi(const char* s1, const char* s2) {
  while (*s1 && (std::toupper(*s1) == std::toupper(*s2))) {
    s1++; s2++;
  }
  return *s1 - *s2;
}

template<typename ElemT>
RandomBC<ElemT>* GetGenClass(const char* gen_name, SeedT seed) {
  if (strcmpi(gen_name, GenNames[Uniform]) == 0) {
    return new RandomDist<ElemT, Uniform>(seed);
  }
  if (strcmpi(gen_name, GenNames[Normal]) == 0) {
    return new RandomDist<ElemT, Normal>(seed);
  }
  MyAbort("Unknown random generator name");  // noreturn!
  return nullptr;   // so never executed.
}
	
  
#endif
