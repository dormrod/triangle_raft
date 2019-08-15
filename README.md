## Triangle Raft

The triangle raft algorithm was originally developed to make hand-built models
of 2D glasses [[1]](https://pdf.sciencedirectassets.com/271563/1-s2.0-S0022309300X04048/1-s2.0-0022309382901065/main.pdf?X-Amz-Security-Token=AgoJb3JpZ2luX2VjEJH%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIBlngSljGlUkTFH54yUE178mACeuABa4uym7N5ddFfyOAiEAsOUoiptcqX6PSGFOY8chuTK%2FpbTegNBpErcYYlqHu7kq2gMIKhACGgwwNTkwMDM1NDY4NjUiDHSw3GrkR2FVigGS2Cq3A5qVBdeIMdsvNk6VI59PcRCpZR1jUxfYxeMNJxZbGrPCoUEFTz1rXHlw0z%2Bj0r6oHxesYkEyclz8kB6NVB%2BH8qeus%2F7othTwkJCOdqa2SsRhXC%2FSCWa53PudrKxFQ4tOQVGE5ZZ%2BohDqHj5lkmk45C0p6gutg7RLAW7r9yskumKw6LjhRXr7j25csj%2BfFUN3qSM9dur8hla0sKdXzKmnnGNgEOq9a2vi6cwUVAdKjES1hCTZ%2BwJQLIVjMwvqBaCUgfMCA24IOEQkvDsgku6lCkGvsbkv8%2B8vULsYLSA2JeVhxlO5ysDAQpIXTkgxsTh9xVnYQw8H%2F83ZGots30A4KgN4060AHSw7ptJzuCkHCIXpnuRo5tPi8JbeEs01dqGg8WL0798qCHXfBrQKvTJ%2BDE%2FuOSThg%2FUzmkDn8cr2RD1i2hw8EH5TSNuBSsxAs0%2BeRuGvfAd5fRphGHkSoitUqEs3gfsmbNDFzXFr3Rv%2Bfv6n6uNIMxG0CX23kAcL6Eub%2FQjSrxNjsqj5Dehb3S6ziM%2BghcckULF45%2BV76wkFyw43yKT9vyrd1kwgR%2F7PWgRqvkLgj2Jehe4wjLTU6gU6tAGUI4WnpVCh41qubTmsBcsYMPoF7KuPfitAzukKGSQO5r9CXQ4pxhUtHmESZC7pxG1pA6mrqwkvZVPjiKQeLk08kXlBPltbu4QI5mIT%2F6XGtnXX%2F2owE7kd8JpwrBANfN3sNHAS%2BZTitOvF%2FjhulslImStYkcMmmHditNyGWMbfhk5rCs6iBitPP55W6QIwkwVbB9Jf1XB%2BgzfSWEc9PKj9UNK5cQ%2BEeXZKTR4LHF1WlA%2F5IIc%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20190815T094349Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY2V6G2QAY%2F20190815%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=efc96c863749e1262fb7182ba3af6192f152fff4fe5b11cd360722b4aca0d165&hash=a978eebb79bc68cd6a6113eefbe42c1a99c544e7fd23a899c3afebd9be95653a&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=0022309382901065&tid=spdf-1eefdc5c-c52e-4013-9cf6-e6c538f703c5&sid=4e47231d31cd204c3759cf153feb0eda1f55gxrqb&type=client).
As part of my PhD I have modified this algorithm into a flexible computational version, to generate networks with tuneable properties [[2]](https://www.tandfonline.com/doi/full/10.1080/00268976.2019.1605099#metrics-content).

The C++ code to generate these structures is included here, along with some python visualisation tools.

For some example outputs please check out the gallery!

![Alt text](./gallery/cover.png "Title")

### Compilation 

Compilation is easiest using CMake.
To compile the code execute the following in terminal:
```commandline
cd src/
cmake .
make
```
As the code is relatively computationally intensive, it is recommended that in the 
```CMakeCache.txt``` you set ```CMAKE_BUILD_TYPE:STRING=Release``` before making.

The generated executable is called ```mx2.x``` (this is for chemistry reasons, the network
in a chemical sense has this formula).

### Input

The parameters for the calculation can be found in ```mx2.inpt``` which are read at runtime.
I would like to say that it is self-explanatory, but as a wise person once said you learn
more about a code if I don't tell you what all the numbers mean...

However, the key parts  are as follows:
```text
1:   /path/to/starting/seed  
```
the algorithm needs a seed to start the calculation, this can be a previous output or
can be generated with the python script```python python_scripts/make_poly_seed.py 6```.
```text:
7:   4 10    ring size limits
```
This pair of numbers determines the size of rings that can be built. 
The lower should be set to 4, the upper can be anything you like >6, but this will slow
the code down!
```text:
8:   2DC/2DS    geometry code
```
If set to 2DC output will be circular, 2DS will be square.
```text:

12:   <temperature>
```
This controls the width of the ring size distribution, the suitable range for which
will depend on the potential model. If this remains unchanged T~10<sup>-5</sup> to 10<sup>-2</sup> should be about right.

### Runtime

Running ```mx2.x``` will read the input file and start the code. A log file, ```mx2.log```, 
will also be produced which is quite verbose. Some things to note:

* The code is not wickedly fast and increasing the ring size limits in particular will slow it down.
* As is the nature of Monte Carlo codes, rarely it will run into irreconcilable cases where the algorithm cannot continue.
In this case the final structure is written. If the temperature is set too high this will happen more often.

### Output and Visualisation

After the code has finished there should be six ```.out``` files, prefixed with your output prefix.
Of these two are probably of the most interest: ```prefix_atoms.out``` gives the atomic number (O or Si), coordination
and coordinates of the final configuration; ```prefix_analysis.out``` contains analysis information
on the ring statistics and bond lengths and angles.

The remaining files specify connectivities. These are needed for visualisation, which can be achieved 
with:
```python
python python_scripts/visualise_network.py prefix <flags>
```
where the choice of flags is as follows:
* -r show rings
* -t show triangles
* -tmx show triangles and atoms
* -s or -S save as .pdf or .png

### Summary

This is a very brief overview of the code, designed to get some output.
Please let me know if you want to use it for research purposes and I will help you 
get to grips with it.

