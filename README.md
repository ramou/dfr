# ThielSort
A competitive implementation of the dfr, written by Larry Thiel. To run, clone this repo:
```
cd ThielSort\dfr
make CannedExtras=July18SimdExtras timing
```
Then run it with 
```
./dfrOpt -n 1000000000 -s -r Uniform
./dfrOpt -n 1000000000 -s -r Normal 1 9223370000000000000 2305840000000000000
./dfrOpt -n 1000000000 -s -r Normal 1 4294970000 1073740000
```
Which runs it for 1 bil vals for Uniform data, a wide normal distribution and a narrow distribution respectively. The "1" says to use a random seed (you could put it after Uniform too), or alternatively you could put a seed for the random generation.

# dfr
Diverting Fast Radix, an incredibly fast algorithm for sorting fixed-digit data.

So, you want to sort things faster because you've realized that under the hood, whatever you're doing is spending way too much time sorting! If your data is even remotely fixed-digit, what we've got here should do the trick. I'll try to make it more user-friendly over time. Sorry that that codebase is absolutely awful right now... but I've been talking a big game for a few years now, and I think I gotta put my money where my mouth is.

## Seeing it run in a simple case
```
git clone https://github.com/ramou/dfr.git
cd dfr
make timing
./perform 1000000
```

This'll give a ton of debug data about where diversion happened and what passes took how long, and how long it took for standard sort to do it. 
If you're super nice, you'll run this on standalone machines for 10<sup>2</sup> till 10<sup>9</sup> and mail me the results, with maybe `/proc/cpuinfo`, `lscpu` or something cool like that. What's neat is I've seen some crazy variation from architecture to architecure that opens the door to the makefile doing a `make install` that squeezes out some crazy extra speed (but only on special cases).

## So, how fast is it?
Check out the [times we've recorded](performance). It's a log scale, so there's a link to the raw data so you can see that.

## How to use it

Until I do a proper `make install` ([#4](https://github.com/ramou/dfr/issues/4)), just make sure you include the [fr.hpp](fr.hpp) file, which contains all the real stuff. 

```C++
template <typename INT, typename ELEM>
void dfr(ELEM *source, auto length)
```

Just call dfr as a templated function with the type INT (which must be the first field) and the type of the overall in-memory object. I currently don't support passing lambdas to get at the fixed-digit key, but [#1](https://github.com/ramou/dfr/issues/1) is about adding that eventually.

If you're sorting an array of `uint64_t`, called `values` that's of length `length`, you would just call
```C++
dfr<uint64_t, uint64_t>(values, length)
```

If you were sorting a bunch of records whose keys were `uint64_t` and whose payload was some other thing and it was of type `ELEM` (with the `uint64_t` key being the first field of `ELEM`) you would call
```C++
dfr<uint64_t, ELEM>(values, length)
```

If things go badly or you'd like me to make it more convenient for you in some way, I'm probably open to making concessions just to get this adopted, so feel free to throw in your two cents.

There are some constants that can yield some improvement via tweaking in the code. [#5](https://github.com/ramou/dfr/issues/5) is about me determining that during `make` and [#4](https://github.com/ramou/dfr/issues/4) is about setting up `make install` so this is even more practical.

## Other folks in the game you should check out
There are two codebases I consider as important to be aware of, both MSD Radix Sorts (I'm in the LSD camp), but both critically relevant. My current codebase is competitive with RADULS2 in the 10s-100s of millions to billions range, so I'll need to squeeze out a few more technical improvements if I want to beat them without non-temporal writes. For smaller inputs (and particularly non-uniform inputs) I start coming out ahead (and will be getting even better). When I say I run faster, I mean A LOT faster. Basically, in the cases where you don't get much for the non-temporal writes, I'm massacaring because LSD Radix Sorts are better, and now I've proven they can divert too (there is an actual proof, I wish I wrote faster/better)! Ska Sort is never competitive with Diverting Fast Radix, but that's not why it's important to recognize. They do a lot of gret stuff that I can learn from. Nobody wants faster unusable code, and studying Ska Sort will make my code more usable.

### Raduls
@marekkokot was friendly, fast and informative in getting back to me and gave me a real competitor to chase. I've specifically avoided non-temporal writes to show that my stuff is competitive without shitting up my code, but he's right to use it, and some day I'll convince a Master's student or someone to add that so I don't need to get my hands dirty... if my dad doesn't get to it first. But hey, I'm also open to collaborating to get it done to appropriately give credit where it is due... I also think you can apply some Fast-Radix-y things to speed up RADULS, but that's just a peripheral thought right now.

RADULS supports multithread, because as an MSD Radix Sort it can. DFR does not. I think the general purpose advantage of this is limited because of the multithread overhead. If you need to sort a million lists, sort each with a single thread and don't worry about unnecessary coordination. If you do need to sort that single list really fast, then multithread becomes way more important.

https://github.com/refresh-bio/RADULS

### Ska Sort
@skarupke has a really neat implementation of a MSD Radix Sort. I think more people should give it the time. My code is easier to use than that RADULS code, but effort has been made to make Ska Sort readily usable in a broad array of cases and I aspire to systematically add the same flexibility to the use of Diverting Fast Radix. A lot of thought has been put into details that I personally haven't cared about, but which I recognize to be really relevant... I hope also to unload much of this on a Master's student as well :D

https://github.com/skarupke/ska_sort

## If you use my code...
Ok, you don't have to. I gotta set the license proper, but I'm pretty down with MIT free. However, if you're using it in a publication, you can reference what I did with a very old version that isn't even half as fast (literally isn't even half as fast):
https://dl.acm.org/citation.cfm?id=2938554
```bibtex
@inproceedings{Thiel:2016:IGL:2938503.2938554,
 author = {Thiel, Stuart and Butler, Greg and Thiel, Larry},
 title = {Improving GraphChi for Large Graph Processing: Fast Radix Sort in Pre-Processing},
 booktitle = {Proceedings of the 20th International Database Engineering \&\#38; Applications Symposium},
 series = {IDEAS '16},
 year = {2016},
 isbn = {978-1-4503-4118-9},
 location = {Montreal, QC, Canada},
 pages = {135--141},
 numpages = {7},
 url = {http://doi.acm.org/10.1145/2938503.2938554},
 doi = {10.1145/2938503.2938554},
 acmid = {2938554},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {algorithm, analytics, big data, graph processing, radix sort},
} 
```

Yes, I'm working with my dad. Yes, it's super cool!


