/* $Header:  $ */ 
#include <limits>

#ifndef RandomToolH
#define RandomToolH

/********************************************************************
 *
 *                            RandomTool.H
 *
 * This namespace provides a selection of random generators that can
 * be selected and configured from command line parameters. For
 * generating random test data.
 ********************************************************************/
#include <bits/stdc++.h>

#include "CmdLineScanner.H"

// RandomTool                                                    RandomTool
namespace RandomTool {

  using EngineT = std::mt19937;
  using SeedT = EngineT::result_type;

  enum Generator { Uniform, Normal };
  const char* GenNames[] { "Uniform", "Normal", "" };

  const char *DefaultGenerator { "Uniform" };
  unsigned long long DefaultSeedSpec { 1 };

  void PrintStringList(const char**const list_v);
  
  // Generic declaration. Only specializations are defined.
  template<typename ElemT, Generator DistId>class RandomDist;

// RandomBC                                                      RandomBC
  // Base class for various random generators. 
  template<typename ElemT> class RandomBC {
    SeedT SeedVal;
    const Generator GenId;
  
  protected:
    RandomBC(Generator rand_dist_id, SeedT new_seed) :
      SeedVal{new_seed}, GenId{rand_dist_id} { }
    
    void BadConfigValues(const char *msg) {
      std::cout << "Error creating configuration for " <<
	"random number generator " <<
	" with (possibly incomplete) configuration :\n";
      PrintConfig();
      std::cout << msg << std::endl;
    }
    void DumpConfig() {
      std::cout << "Random number generator config: Generator = " <<
	GenNames[GenId] << " Seed = " << SeedVal << ' ';
    }
public:
    SeedT Seed() { return SeedVal; }  // Return the current seed.
    void Seed(unsigned new_seed) { SeedVal = new_seed; }
    // Printable name of distribution
    const char *WhoAmI() { return GenNames[GenId]; }

    // Fill data[] with n random values. 
    virtual void MakeRandomVec(ElemT *data, long n) = 0;  
    virtual void Configure(CmdLineScanner &cls) = 0;
    virtual void PrintConfig() = 0;
  
  };  // RandomBC

// Uniform                                                        Uniform
  template<typename ElemT>class RandomDist<ElemT, Uniform> :
    public RandomBC<ElemT> {
    static const Generator DistId = Uniform;
    using BC = RandomBC<ElemT>;
    ElemT Min = 0;
    ElemT Max = std::numeric_limits<ElemT>::max();

  public:

    RandomDist(SeedT seed_val) : BC(DistId, seed_val) {}
    ~RandomDist() {}

    void MakeRandomVec(ElemT *data, long N) {
      std::mt19937 generator(this->Seed());
      std::uniform_int_distribution<ElemT> distribution(this->Min, this->Max);
      for(unsigned long long i = 0; i<N;i++)data[i]=distribution(generator);
    }
    void Configure(CmdLineScanner &cls) {
      if (!cls.ExpectUnsignedValue(&Min)) return;  // Default config
      if (!cls.ExpectUnsignedValue(&Max))  // Can't give just the first!
	MyAbort("Invalid config for uniform dist");
    }
    static void ConfigHelp() {
      const char* HelpText
	[] { "Configuration for Uniform dist: 2 values of type ElemT:\n  ",
	  "MinValue, MaxValue. Defaults: 0, ElemT::max\n", ""};
      PrintStringList(HelpText);
    }
    void PrintConfig() {
      this->DumpConfig();
      std::cout << "Min = " << Min << " Max = " << Max << std::endl;
    }
			    
  };
// Normal                                                           Normal
  template<typename ElemT>class RandomDist<ElemT, Normal> :
    public RandomBC<ElemT> {
    static const long long DefaultAvg { 1l<<32 };
    static const long long DefaultDev { 1<<30 };

    static const Generator DistId = Normal;
    using BC = RandomBC<ElemT>;
    long long Avg=DefaultAvg;
    long long Deviation=DefaultDev;

  public:
    
    RandomDist(SeedT seed_val) : BC(DistId, seed_val) {}
    ~RandomDist() {}
    
    void MakeRandomVec(ElemT *data, long N) {
      std::mt19937 generator(this->Seed());
      std::normal_distribution<> distribution(this->Avg, this->Deviation);
      const double bottom = 1.0;
      const double top = std::numeric_limits<ElemT>::max();
      for(unsigned long long i = 0; i<N;i++) {
	do {
	  double t = distribution(generator);
	  if ((t < bottom) || (t > top)) continue;
	  data[i]= t;
	  break;
	} while (true);
      } 
    }
    void Configure(CmdLineScanner &cls) {
      if (!cls.ExpectUnsignedValue(&Avg)) return;  // Default both
      if (!cls.ExpectUnsignedValue(&Deviation)) return;  // Default deviation
    }
    static void ConfigHelp() {
      const char* HelpText {
	"Configuration for Normal dist: 2 values of type double:\n  "
	  "AverageValue, Deviation. Defaults: " } ;
      std::cout << HelpText << DefaultAvg << ", " << DefaultDev << "\n";
    }
    void PrintConfig() {
      this->DumpConfig();
      std::cout << "Avg = " << Avg << " Deviation = " <<
	Deviation << std::endl;
    }
  };

  void PrintDistHelps() {
    RandomDist<unsigned, Uniform>::ConfigHelp();
    RandomDist<unsigned, Normal>::ConfigHelp();
  }

#include "RandomToolHelper.H"

  /**
   * Extract the random dist spec from the command line arguments. The
   * start of the random spec is at arg_ndx. On return arg_ndx has
   * been incremented to point to the first parameter after the random
   * spec, or ie == argc if there is none. Return a pointer to a new
   * object of the requested generator class, which has been
   * configured according to the command line configuration.
   * Abort with a message if an invalid spec is found.
   **/
  template<typename ElemT>
  RandomBC<ElemT>* GetGenClass(CmdLineScanner *cls) {
    if (!cls) {  // Use defaults.
      return GetGenClass<ElemT> (DefaultGenerator, GetSeed(DefaultSeedSpec));
    }
    const char *genName{nullptr};
    if (!cls->ExpectString(&genName)) MyAbort("Bad/missing gen name");
    unsigned long long t = 0;
    if (!cls->ExpectUnsignedValue(&t)) {  // Default seed, default config.
      return GetGenClass<ElemT>(genName, GetSeed(DefaultSeedSpec));
    }
    RandomBC<ElemT>* rbc = GetGenClass<ElemT>(genName, GetSeed(t));
    if (!cls->IsEnd()) {  // Might be a configuration.
      rbc->Configure(*cls);
    }
    return rbc;
  }

};  // RandomTool


#endif
