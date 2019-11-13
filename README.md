# dfr
Diverting Fast Radix, an incredibly fast algorithm for sorting fixed-digit data. So, you want to sort things faster because you've realized that under the hood, 
whatever you're doing is spending way too much time sorting! If your data is even remotely fixed-digit, what we've got here should do the trick. I'll try to make 
it more user-friendly over time. Sorry that that codebase is absolutely awful right now... but I've been talking a big game for a few years now, and I think I 
gotta put my money where my mouth is.
## Seeing it run in a simple case
``` git clone https://github.com/ramou/dfr.git make perform ./perform 1000000 ``` This'll give a ton of debug data about where diversion happened and what passes 
took how long, and how long it took for standard sort to do it. If you're super nice, you'll run this on standalone machines for 10<sup>2</sup> till 
10<sup>9</sup> and mail me the results, with maybe `/proc/cpuinfo` or something cool like that. What's neat is I've seen some crazy variation from architecture to 
architecure that opens the door to the makefile doing a `make install` that squeezes out some crazy extra speed (but only on special cases).
## How to use it
Until I do a proper `make install` ([#4](https://github.com/ramou/dfr/issues/4)), just make sure you include the [fr.hpp](fr.hpp) file, which contains all the 
real stuff. ```C++ template <typename INT, typename ELEM> void dfr(ELEM *source, auto length) ``` Just call dfr as a templated function with the type INT (which 
must be the first field) and the type of the overall in-memory object. I currently don't support passing lambdas to get at the fixed-digit key, but 
[#1](https://github.com/ramou/dfr/issues/1) is about adding that eventually. If you're sorting an array of `uint64_t`, called `values` that's of length `length`, 
you would just call ```C++ dfr<uint64_t, uint64_t>(values, length) ``` If you were sorting a bunch of records whose keys were `uint64_t` and whose payload was 
some other thing and it was of type `ELEM` (with the `uint64_t` key being the first field of `ELEM`) you would call ```C++ dfr<uint64_t, ELEM>(values, length) ``` 
If things go badly or you'd like me to make it more convenient for you in some way, I'm probably open to making concessions just to get this adopted, so feel free 
to throw in your two cents. There are some constants that can yield some improvement via tweaking in the code. [#5](https://github.com/ramou/dfr/issues/5) is 
about me determining that during `make` and [#4](https://github.com/ramou/dfr/issues/4) is about setting up `make install` so this is even more practical.
## Other folks in the game you should check out
There are two codebases I consider as important to be aware of, both MSD Radix Sorts (I'm in the LSD camp), but both critically relevant. My current codebase is 
competitive with RADULS2 in the hundreds of millions to billions range, so I'll need to squeeze out a few more technical improvements if I want to beat them 
without non-temporal writes. For smaller inputs (and particularly non-uniform inputs) I start coming out ahead (and will be getting even better). When I say I run 
faster, I mean A LOT faster. Basically, in the cases where you don't get much for the non-temporal writes, I'm massacaring because LSD Radix Sorts are better, and 
now I've proven they can divert too (there is an actual proof, I wish I wrote faster/better)! Ska Sort is never competitive with Diverting Fast Radix, but that's 
not why it's important to recognize. They do a lot of gret stuff that I can learn from. Nobody wants faster unusable code, and studying Ska Sort will make my code 
more usable. @marekkokot was friendly, fast and informative in getting back to me and gave me a real competitor to chase. I've specifically avoided non-temporal 
writes to show that my stuff is competitive without shitting up my code, but he's right to use it, and some day I'll convince a Master's student or someone to add 
that so I don't need to get my hands dirty... if my dad doesn't get to it first. But hey, I'm also open to collaborating to get it done to appropriately give 
credit where it is due... I also think you can apply some Fast-Radix-y things to speed up Raduls, but that's just a peripheral thought right now. 
https://github.com/refresh-bio/RADULS @skarupke has a really neat implementation of a MSD Radix Sort. I think more people should give it the time. My code is 
easier to use than that RADULS code, but effort has been made to make Ska Sort readily usable in a broad array of cases and I aspire to systematically add the 
same flexibility to the use of Diverting Fast Radix. A lot of thought has been put into details that I personally haven't cared about, but which I recognize to be 
really relevant... I hope also to unload much of this on a Master's student as well :D
https://github.com/skarupke/ska_sort
