# Mandelbrot-set-topology
Mandelbrot set topology

How to read informations from images ?

How to describe hyperbolic component of Mandelbrot set ? ( escpecially [island](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/island_t) )
* period of main pseudocardioid
* center of main pseudocardioid
* window ( radius and center) of parameter plane
* angles of external rays that land on cusp of the pseudocardiod
* size of main pseudocardioid
* orientation
* [distortion](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Mandelbrot_set/mset_distortion)

# To do

For each island: 
* find period 
* find center of component ( nucleus) = pseudocardioid_nucleus
* find shape ( psudocircle or pseudocardioid)
* find cusp of pseudocardioid = pseudocardioid_cusp
* find root point between 2 main componnents = pseudocardioid_root_half
* find center of second component ( nucleus) = second_bulb_nucleus
* compute [distortion of the island](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Mandelbrot_set/mset_distortion)



# Example islands

Here are few examples of [islands](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/island_t). 
* first image ( LastIteration = 108):  shows interior points
* second image (Period = 107): Only main (pseudo)cardioid of period p  and main component of period 2p is drawn 
* for more description see the [output the new txt file m.txt](./src/cli/new/m.txt) and [output the old txt file m.txt](./src/cli/old/m.txt)

```c
#define kMax 21 // number of examples, see plane_examples

//plane : plane_center_x,  		plane_center_y,  		plane_radius, 		period

double plane_examples[kMax][4] = {
	{-0.4,				+0.0,				0.8,			1}, // k = 0
	{+0.29254,			-0.01497, 			0.00015,		16}, 
	{-1.763,  			+0.0,				0.016,			3}, 
	{-0.15842, 			+1.03335, 			0.008,			4},  
	{+0.358431,			+0.643507,			0.005,			5},  
	{+0.442990,			+0.374,				0.003,			6},  
	{+0.432259,			+0.2275,			0.002,			7}, 
	{+0.404879,			+0.146216,			0.0015,			8}, 
	{+0.378631,			+0.098841,			0.001,			9}, 
	{+0.356854, 			+0.069659,			0.0007,			10},
	{+0.339454,			+0.050823,			0.0005,			11}, // k = 10
	{+0.325631,			+0.038164,			0.0005,			12},
	{+0.260270,		 	+0.00167,			0.00002, 		32 },
	{+0.2524945,			+0.0001973,			0.0000025,		64},
	{+0.25061329,			+0.00002399,			0.0000003,		128},
	{+0.250151979,			+0.000002959,			0.000000036,		256},
	{+0.250037823,			+0.0000003673,			0.000000004,		512}, 
	{+0.2500094342,			+0.000000045694786520646,	0.0000000005,		1024}, //
	{+0.25000235583,		+5.701985912706845832e-09,	0.00000000007,		2048}, // 
	{+0.2500005886144,		+0.0000000007122,		0.0000000000083,	4096}, //
	{+2.500001471109009610e-01,	+8.897814201389663379e-11,	1.1e-12,		8192}  //  long time 
	 
};

```
## Period 1
peroid = 1 = whole Mandelbrot set  ( using LastIterarion method with interior detection)

![](./png/0_period_1_LastIteration.png "period = 1 = whole Mandelbrot set ") 

Here unknown pixels are marked by red big pixels. They are boundary points. 

![](./png/0_period_1_LastIteration_unknown.png "period = 1 = whole Mandelbrot set with unknown pixels marked red ") 

![](./png/0_107.png "description") 


With important points: pseudocardioid_nucleus, pseudocardioid_cusp, pseudocardioid_root_half, second_bulb_nucleus  

![](./png/0_114.png "description") 


## Period 16 island

![](./png/1_period_32_LastIteration.png "description") 

Here are 2 main components : period 16 pseudocardioid and period 2*16 pseudocircle. Note that the image differs from that computed with LastIteration method ( above). 
The components are smaller and do not touch in root points

![](./png/1_107.png "description") 

Here is the comparison using ImageMagic ( [see bash file](./src/cli/old/m.sh))

![](./png/diff1.png "description") 

Here is a comparison between boundaries of LastIterarion method and Period method. The difference is big. It was caused by to low value of iMax = 2000; ( see function GivePeriodByIteration line 676 )

![](./png/1_111.png "description") 

After increasing iMax (10 000 and 50 000) ( and orbit[OrbitLength]; // length(orbit) = iMax + 1 ) new image looks better 

![](./png/1_111_10000.png "description") 

![](./png/1_111_50000.png "description") 


With 
* important points: pseudocardioid_nucleus, pseudocardioid_cusp, pseudocardioid_root_half, second_bulb_nucleus  
* 2 axis, one thru pseudocardioid_nucleus and pseudocardioid_cusp, second thru  pseudocardioid_nucleus and pseudocardioid_root_half
* distorsion angle is denoted with arc

![](./png/1_114.png "description") 


## Period 3

Period 3 island looks like whole Mandelbrot set  

![](./png/2_period_3_LastIteration.png "description") 

![](./png/2_107.png "description") 

With important points, axes and distortion angle 

![](./png/2_114.png "description") 

## Period 4

![](./png/3_114.png "description") 

## Period 5

![](./png/4_114.png "description") 

## Period 6


![](./png/5_114.png "description") 

## Period 7

![](./png/6_114.png "description") 

## Period 8

![](./png/7_114.png "description") 

## Period 9


![](./png/8_114.png "description") 

## Period 10


![](./png/9_114.png "description") 

## Period 11


![](./png/10_114.png "description") 

## Period 12 

![](./png/11_114.png "description") 


## Period 32 

![](./png/12_114.png "description") 


## Period 64 

![](./png/13_114.png "description") 


## Period 128

![](./png/14_114.png "description") 


## Period 256

![](./png/15_114.png "description") 


## Period 512 

![](./png/16_114.png "description") 


## Period 1024 

![](./png/17_114.png "description") 

## Period 2048 

![](./png/18_114.png "description") 


## Period 4096 

![](./png/19_114.png "description") 


## Period 8192 

![](./png/20_114.png "description") 

# Files
* [m.c](./src/cli/new/m.c) - one file c program
* [m.sh](./src/cli/new/m.sh) - bash file used to compile, run the m.c and image manipulations using ImageMagic console programs
* [Makefile](./src/cli/new/Makefile) file for make which runs all
* [output txt file m.txt](./src/cli/new/m.txt) - to read about the results
* [old code](./src/cli/old/) 
* [images are in the png directory](./png/) 


# Algorithms
* [atom domains for period  (= period domains) of Mandelbrot set hyperbolic components](https://commons.wikimedia.org/wiki/File:Mandelbrot_Atom_Domains_Animation.gif)
* period of Mandelbrot set hyperbolic componnets
* [multiplier map](https://commons.wikimedia.org/wiki/File:Mandelbrot_set_-_multiplier_map.png)
* [interior detection](https://commons.wikimedia.org/wiki/File:Mandelbrot_set_with_Interior_detection_method.png)
* [The Lyapunov exponent](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Mandelbrot_set_interior)
* [Interior distance estimation](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/demm#Interior_distance_estimation) - DEM 
* [rendering-mandelbrot-set-mu-units](https://fractalforums.org/code-snippets-fragments/74/rendering-mandelbrot-set-mu-units/3485): How to extract a mu-unit of a given period using distance estimation colouring, pruning off the outer filaments?
* [mu-atom ](https://mathr.co.uk/mandelbrot/mu-atom/) - mu-atom mapping: period p hyperbolic components of the Mandelbrot set can each be mapped conformally to the unit disc, by the derivative d/dz of the periodic limit cycle where f_c^p(z_0) = z_0.
* [period-doubling-in-minibrots](https://fractalforums.org/noobs-corner/76/period-doubling-in-minibrots/3990)
* [shape estimation]()


Description
* all calculations are in a point-sampling, numerical ( double precision), non-rounding controlled manner
* execution time of DrawImage for LastIteration = 31 	for Period = 525 	MultiplierMap = 5 seconds for example 2 ( period 3)





# Key words
* [distortion](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Mandelbrot_set/mset_distortion)
* [computational topology](https://en.wikipedia.org/wiki/Computational_topology)
* [digital topology](https://en.wikipedia.org/wiki/Digital_topology)
* 2D 
* [raster graphic](https://en.wikipedia.org/wiki/Raster_graphics)
* computer graphic
* [binary image](https://en.wikipedia.org/wiki/Binary_image)
* [Connected-component labeling](https://en.wikipedia.org/wiki/Connected-component_labeling)



# Git

create a new repository on the command line
```
echo "# " >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:adammaj1/Mandelbrot-set-topology.git
git push -u origin main
```


## Repo

Change:
* in general settings
  * add Social Preview Image ( Images should be at least 640??320px (1280??640px for best display))
* in repository details ( near About) add
  * description
  * website 
  * Topics (separate with spaces) 
  

Local repository

```
~/Dokumenty/Mandelbrot-set-topology/ 

```





## Subdirectory

```git
mkdir png
git add *.png
git mv  *.png ./png
git commit -m "move"
git push -u origin main
```
then link the images:

```txt
![](./png/n.png "description") 

```

to overwrite

```
git mv -f 
```

```
git mv ./src/*.c ./src/modified/bash/
git mv ./src/*.sh ./src/modified/bash/
gitt mv ./src/Makefile ./src/modified/bash/
```



## Github
* [GitHub Flavored Markdown Spec](https://github.github.com/gfm/)
* [md cheat sheet](http://mdcheatsheet.com/)
* [CommonMark Spec](https://spec.commonmark.org)
* [Markdown parser ](https://markdown-it.github.io/)

