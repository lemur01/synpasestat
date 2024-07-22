# synpasestat

## Analysis of GABA/glutamate signal.

The script SynapseStats performs the analysis of images of neurons grown of substrates of different stifness.
Input: folder containing lif files with four colour experiments: neurofilaments (channel 1), glutamate (channel 2), neuroligin (channel 3), and GABA transporter (VGAT, channel 4).
Output: result folder with counts for the above signals (total, close to one/two/all other colours respectively) and length of detected process.

 ## Parameters:
- thDist    (default 2)       - distance threshold in pixels (for signal to be considered close)
- thRemove  (default 150)     - process fragments smaller than the threshold are removed
  
<img src=src/Synapse.png width="300" height="300">

### External dependency
bfmatlab
https://www.openmicroscopy.org/bio-formats/

Melissa Linkert, Curtis T. Rueden, Chris Allan, Jean-Marie Burel, Will Moore, Andrew Patterson, Brian Loranger, Josh Moore, Carlos Neves, Donald MacDonald, Aleksandra Tarkowska, Caitlin Sticco, Emma Hill, Mike Rossner, Kevin W. Eliceiri, and Jason R. Swedlow (2010) Metadata matters: access to image data in the real world. The Journal of Cell Biology 189(5), 777-782. doi: 10.1083/jcb.201004104



