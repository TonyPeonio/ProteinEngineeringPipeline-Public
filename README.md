# Customized End-to-End Pipeline for Protein Engineering Beginners

**Hello!** Whether you're a seasoned developer, a biology student, or just someone curious about computational biology, welcome. If you have any questions, run into issues, or want to suggest improvements, please feel free to open an issue or reach out. I'd love to connect and discuss the pipeline!

An automated, start-to-finish pipeline for designing novel proteins, built specifically for people who want to explore protein engineering without needing experience in biology (though a willingness to learn is important!).

**The Philosophy:** I built this as an undergraduate Computer Science student who wanted to learn about protein engineering but struggled with the sheer amount of domain knowledge usually required. This pipeline is designed to maximize reliance on compute power and minimize reliance on human input. The tedious, knowledge-heavy steps (like identifying hotspots) are not completely automated, but have custom methods to attempt to decrease the necessity for field expertise. 

My goal is simple: make protein engineering accessible to anyone with a GPU (or anyone willing to spin up Google Cloud) and an interest in the field. **This pipeline (and the attached resources) were enough to win a University Protein Engineering Competition (open to both undergraduates and graduates) with prize money.**

I know that READMEs are easy to skip, but I strongly encourage that you read through this in its entirety. I did my best to keep it brief while covering everything that you will need to know to use this pipeline.

**Disclaimer**: Though I have learned a lot throughout the creation of this pipeline, there are components of this readme/pipeline that might be incorrect in some regards. While I have done my utmost to prevent this and discuss the code with others, there may still be mistakes.

---

## 1. Setup and Installation

This pipeline brings together several powerful tools (RFdiffusion, ProteinMPNN, LocalColabFold, and PyRosetta). Because machine learning environments can be tricky, the setup is split into two parts: an automated script for the basics, and manual instructions for the hardware-specific libraries. 

Despite trying, fully automating the C++ and CUDA bindings across all hardware was beyond my current skillset, so some manual work is left for you in Step 2. However, this README (and the official RFdiffusion README) are intended to be strong resources for navigating this part. The rest of the setup is completely automated for you!

### Step 1: Get the Code and Run the Setup Script
First, download the repository to your machine and run the automated setup script. This script will download the required repositories, set up basic conda environments, and download the necessary model weights.

```bash
# Clone the repository
git clone https://github.com/TonyPeonio/ProteinDesignChallenge.git

# Navigate to the scripts folder, make it executable, and run
cd ProteinDesignChallenge/scripts_and_helpers
chmod +x setup.sh
./setup.sh
```

Step 2: Hardware-Specific Configuration (RFdiffusion)
Because PyTorch and DGL require specific C++ binaries depending on your graphics card, so you must install them manually. Run the following commands based on your system. (Note: If these specific versions do not work for your hardware, please refer to the official RFdiffusion README for troubleshooting. If you are using a new version of cuda (necessary with GPUs like the rtx 5070 oc), it might be necessary to modify the miniconda source code by going into /miniforge3/envs/env_rfdiffusion/lib/python3.9/site-packages/dgl/graphbolt/__init__.py and commenting out the bottom line (load_graphbolt()), which has the potential to get rid of one of the errors you may encounter while trying to balance new pytorch with old things like dgl and graphbolt).

1. Activate the environment:

```bash
conda activate env_rfdiffusion
```
2. Install PyTorch and DGL (Choose ONE based on your hardware):

For Modern NVIDIA GPUs (RTX 30-series, 40-series, 50-series / CUDA 12.1+ (possibility for new GPUs to require nightly torch)):

```bash
pip install torch==2.8.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128
pip install dgl -f https://data.dgl.ai/wheels/cu121/repo.html
```

This might not be necessary on all GPUs, but the only way I could get my RTX 5070 OC working was by going into ~/miniforge3/envs/env_rfdiffusion/lib/python3.9/site-packages/dgl/graphbolt/__init__.py and modifying the final line, as noted in the instructions above.

For Legacy NVIDIA GPUs (RTX 20-series or older / CUDA 11.8):
```bash
pip install torch==2.4.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install dgl -f https://data.dgl.ai/wheels/cu118/repo.html
```
For CPU-only (Mac or no GPU):
```bash
pip install torch==2.4.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install dgl -f https://data.dgl.ai/wheels/repo.html
```
3. Install core dependencies:

```bash
pip install decorator scipy hydra-core omegaconf pyrsistent "torchdata<0.10.0" pandas
```
4. Compile the SE3Transformer:
```bash
cd ../RFdiffusion/env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../..
pip install -e .
```

## 2. Pre-Run Intuition and Considerations (Before You Code)

Before you fire up the pipeline, there are a few practical things you should consider about the protein you want to target. You don't need a biology degree, but keeping these constraints in mind will save you time and compute power. If you would like more information, there are many excellent resources online, some of which are linked in scripts_and_helpers/docs/useful_info.txt. All configs discussed are set in scripts_and_helpers/scripts/run_all.sh:

* **VRAM and Protein Size:** AI models are extremely hungry for VRAM. If you try to model a massive complex (for example, I initially tried to model the entire Nitrogenase complex), your GPU will run out of memory and crash. Even if it doesn't crash and overflows to standard RAM, it is **much** slower (something to watch out for). **The Fix:** Find your target protein on the RCSB PDB, look at its 3D structure, and identify just the specific chain(s) you actually want to bind to. Crop the rest out. PyMOL is an excellent resource for viewing and manipulating proteins, and one that I would strongly recommend for every single step in this pipeline (except ProteinMPNN). 
* **Inhibitory vs. Non-Inhibitory:** What is your goal? 
    * **Non-Inhibitory**: If you just want to attach a biological "tag" to the protein without breaking it, you can let the AI find any stable binding pocket on the surface, though you do need to have a few considerations:
        * Ensure that it does not prevent a docking protein that brings resources.
        * Ensure that it does not warp the shape of the protein to the point that it cannot function.
        * Ensure that it does not block the active sites of the protein
    * **Inhibitory**: If you want it to be inhibitory, just do one of the things you try to prevent with a non-inhibitory binder.
* If the protein is small enough that it fits within the VRAM and you are a beginner, do not worry about cropping the protein. Go and set the SERINE_PAINTING config to 0. Set the CORE_PROTECTION to 1. 
* If the protein is large, choose the chain (or chains, depending on chain size and VRAM) and enter the chains and amino acid ranges into the AA_RANGE config. If there are multiple chains or segments, add a /0 between each one.

---

## 3. How to Run the Pipeline

The beauty of this pipeline is that the entire automated pipeline is controlled from a single file: `run_all.sh`.

You do not need to run a dozen different scripts. You simply open `run_all.sh`, fill out the Master Configuration section at the top, and let it run. The script will automatically fetch your target PDB, clean it, generate binder shapes (RFdiffusion), assign sequences (ProteinMPNN), and validate whether they actually fold (LocalColabFold). After this part, understanding some biology is helpful, but I will walk you through the PyRosetta script later in the README. Feel free to get the run_all.sh started after filling out the config before you begin reading the PyRosetta scripts.

**To start the process:**
1. Open `ProteinDesignChallenge/scripts_and_helpers/scripts/run_all.sh` in a text editor.
2. Edit the variables in the `MASTER CONFIGURATION` section (explained below).
3. Save the file and run it:
```bash
./run_all.sh
```

## 4. Understanding the Configuration Settings

There are a wide variety of configuration settings that mean incredible customization opportunities, all of which are discussed in detail below. To get a general intuitive idea for the pipeline (which is necessary to understanding the configurations), I will give a brief conceptual overview of how the pipeline works. I will also give time estimations for each step as a frame of reference, but this could vary wildly based on hardware and problem.

First, RFdiffusion creates blank pdbs. It directly or indirectly utilizes the following configs:
**PDB_ID,DESIGNABLE_CHAINS,BINDER_LENGTH,NUM_DESIGNS,AA_RANGE,HOTSPOTS**.
RFdiffusion uses a noise-based AI model to attempt to create novel binders. The functionality of RFdiffusion is to create backbones for the binder that you want to create. It is essentially just a shape, in PDB format, with x,y,z locations that fit the parameters that you gave it, but it does not have any real idea of how it is created. For example, there is no valid amino acid sequence. On my hardware, on average, one backbone took about 4 minutes to generate for me.

This is where ProteinMPNN comes in. ProteinMPNN takes in the blank PDBs that RFdiffusion output, and turns them into a probable fasta sequence (a text sequence of amino acids, using one letter acronyms). On average, one fasta sequence took about 12 seconds to generate for me.

To confirm that this fasta sequence actually binds as intended, it needs to be run through a folding software. LocalColabFold is derived from AlphaFold2, and does an excellent job of locally analyzing your fasta sequences and turning it back into a pdb. This can output a number of important metrics (iptm, ptm, plddt) that will be discussed in more detail later. On average, this took about 5 minutes per fasta sequence with 5 models.

Finally, after all of these steps, there are provided tools to manually manipulate the proteins using PyRosetta, which has the power to do chemical interface analysis, mutagenesis, and relaxations, all of which are incredibly useful tools. Using these tools will be discussed in more detail later. This step can wildly vary in time requirements, depending on the depth of analysis, ranging from 10 minutes to 8 hours.

Once you have your final results, it is recommended to enter the FASTA sequences produced by the PyRosetta scripts to confirm their success into AlphaFold3. This will again produce useful metrics (iptm, ptm, plddt). This step can take between 0.5 and 12 hours, depending on server availability.

I will also give examples of the configuration that I used to complete the Protein Engineering Competition mentioned above. The goal of the competition was to create a non-inhibitory binder for nitrogenase.

### The Configuration Parameters

* **PDB_ID**: Find your protein in the RCSB Protein Data Bank. Enter the 4-character code. The pipeline will automatically handle importing and cleaning the PDB. For example, my configuration was 7ut8.

* **DESIGNABLE_CHAINS**: Only plug in chains that you want to be modified. Do not enter all of the chains. For example, my configuration was B (I cropped only chain A out of the nitrogenase structure, meaning B was defined to be my binder). 

* **BINDER_LENGTH**: The length of the binder that you are trying to create. For example, my configuration was 55-65 (after lots of iteration and learning about what was best for my specific problem). 

* **NUM_DESIGNS**: The number of iterations that you want to run. I would recommend doing batches of 200-300 until you find a good configuration set, unless you have a very high amount of computer power. For example, my configuration on my first runs were only on the order of 100, while I ended on the order of 1000 when I found successful configs.

* **AA_RANGE**: The chains and amino acids you would like to pull from the pdb to actually use in the simulations. As mentioned above, this pipeline does allow for cropping. For example, my configuration was "A5-480", which is the entire chain A. Something important to note is that chains do not always start at 0 or 1, and occasionally, there are situations where some numbers are skipped. If this occurs and dashes appear in the fasta sequence (essentially a text sequence with a single letter representing each amino acid), the script automatically replaces them with G, as it is highly flexible. If this is not the preferred behavior, modify 3-2_mpnn_sequence_filtering.py. I have not encountered this often, but wanted to give a warning.

* **HOTSPOTS**: This will be discussed in detail later. This pipeline has a special model that I came up with (might also exist elsewhere, but I independently created this method), which drastically reduces the reliance on the necessity for expertise. This, along with the BINDER_LENGTH and AA_RANGE, are likely the most important configuration settings along with the most difficult to derive. How to find all of them will be discussed in detail below. For example, my configuration was "A30,A32,A34,A36". 

* **SEQ_PER_BACKBONE**: The number of fasta sequences that proteinMPNN creates for each RFdiffusion backbone. Raising this number could increase the structural stability of the binder but decrease the interface quality. For example, I commonly just used 1, because I still got high structural stability, and the main thing I was trying to maximize was interface quality.

* **MPNN_TEMP**: The amount of creativity ProteinMPNN uses when creating fasta sequences. Higher numbers are more creative. For example, I commonly used 0.1, though I did experiment with up to 0.3.

* **WILDCARDS**: Not something to worry about if **SEQ_PER_BACKBONE** is 1. If so, leave **WILDCARDS** at 0. If not, how the script works is to get the best sequence for each of the backbones, then get an extra X wildcards. The reasoning behind the method of not just taking the top X wildcards and ignoring if one is the best for a given sequence is to maximize diversity. For example, I commonly had this value at 0.

* **SERINE_PAINTING**: This was one of the highly custom parameters added to allow for increased ability to crop without inaccurate results. This parameter is highly related to **CORE_PROTECTION**, so ensure that you read both descriptions. If there is no cropping, set this to 0. If there is cropping, continue reading. When a protein is cropped, some parts that were only exposed to other chains, and often have extremely attractive energy pockets, are now exposed. However, a protein getting stuck in this energy well is unfair, as it would not be available in the lab to bind to, as it would be covered by other chains. Therefore, the idea is to paint over it with Serine, which is an unattractive amino acid. This will search for all atoms within 8 angstroms of where a ghost chain (ghost chain being a chain that was cropped away), and mutate it to Serine. If **SERINE_PAINTING**=1, all of those atoms will become Serine (with the exception of **CORE_PROTECTION**, explained below).  If **SERINE_PAINTING**=0.5, only 50% of those atoms become **Serine**, determined by sorting the list by importance. The goal is to find a **SERINE_PAINTING** value as high as possible while not disrupting the structure of the protein to maximize the number of binders that bind elsewhere and prevent false negatives. I invented these parameters after the competition ended, so I cannot provide useful values. While I was using Serine Painting, I was using hardcoded values that do not easily translate to these configs.

* **CORE_PROTECTION**: As you might expect, painting a surface with Serine could lead to problems with the structure of the protein. If too much becomes Serine, it loses its structure, warps, and begins creating false positives/negatives with the binders. A good way to determine if this is happening is to align the chain with the wildtype chain in PyMOL (there are good online resources for this), and the RMSD will be printed in the terminal. For a decently sized chain, less than 1.0 means it is fairly safe, though the lower the better. If it is above 1, become concerned that the protein is being warped too much. This is where **CORE_PROTECTION** comes in. Assume **CORE_PROTECTION** is set to 0.7. Then, all of the amino acids are ranked by how structurally important they are, and the most important 70% are now allowed to be turned into Serine. I invented these parameters after the competition ended, so I cannot provide useful values. While I was using Serine Painting, I was using hardcoded values that do not easily translate to these configs.

* **NUM_RECYCLES**: I strongly recommend keeping this value at 12 unless you have a good reason to change it. This is the maximum number of iterations for which LocalColabFold can try to wiggle around the protein to see if it can get a better fit. While this may seem like a large number, the script is set to automatically terminate if it does not improve by the **RECYCLE_TOLERANCE** from one iteration to the next, so it very rarely reaches 12. I exclusively kept it at 12 once I learned how it actually worked and how to implement tolerance, and my results were drastically better than when I started with it at 2.

* **RECYCLE_TOLERANCE**: This is the amount by which LocalColabFold needs to modify the fold between consecutive iterations to try again. 0.5 is an excellent value and I would recommend leaving it here unless you have cause to change it.

* **NUM_MODELS**: LocalColabFold comes with 5 different models that can run on your fasta sequence. Some are more or less pessimistic about binding locations, folding structures, and more. This diversity is extremely powerful for helping to recognize when a binder is useful. It will not always be useful on all models, so I would recommend keeping this at 5 if possible. It does take 5x as long as 1 model, but I would still recommend keeping it at 5, especially because 5 models on LocalColabFold only takes about as long as one backbone on RFDiffusion (in my experience). 

**Congratulations! You are ready to begin your first run!** If you are confident about the config, fill out the run_all.sh with the appropriate config for your needs, navigate to the scripts_and_helpers/scripts directory, and then run the script with 
```bash
./run_all.sh
```

### Configuration Parameters - Extra Help
If you are not confident in what to put for the config, I listed a process below for finding a good HOTSPOT config. Most other configs you can figure out just based on your problem, repeated experimentation, or the examples I provided for my work in the Configuration Settings Section.

* If you are **CROPPING**:
  * Install PyMOL for the following steps (or another software that can calculate RMSD when aligning two chains)
  * Test out different values of SERINE_PAINTING and CORE_PROTECTION until you can find a value that is around 0.8 RMSD when aligning the chains in PyMOL
  * Run 500 iterations with a random valid hotspot (the hotspot is valid if it is an existing amino acid in the chains passed into AA_RANGE). Only choose one amino acid as a hotspot, not multiple
  * Open scripts_and_helpers/1_colabfold_sequence_filtering.ipynb
    * Do not run the first cell
    * Run the second cell and wait for the calculation to finish
      * When finished, it will show a histogram that is a visual representation of where the binder wants to bind
        * The reasoning behind this histogram is working with the protein instead of against the protein
        * Using this knowledge, examine the high energy points on the graph. This is where the protein wants to bind
        * If you are seeking an inhibitory protein, find a spot that blocks some sort of active site. This varies widely based on the problem so I cannot give more details
        * If you are seeking a non-inhibitory (more difficult), consider all possible inhibitory actions and find the attractive spot that affects the protein the least. This is a great spot to refer to external resources like youtube/papers to find out what makes a spot inhibitory/non-inhibitory.
    * Using the outputs from the Jupyter notebook, update the hotspots to the site that you determined was best for your needs
    * I would recommend creating a wide range of BINDER_LENGTH so the RFdiffusion has freedom to explore until you learn what is best
    * Run again and see how it goes!
* If you are **NOT CROPPING**:
  * Set SERINE_PAINTING and CORE_PROTECTION to 0 and 1, respectively
  * Run 500 iterations with a random valid hotspot (the hotspot is valid if it is an existing amino acid in the chains passed into AA_RANGE). Only choose one amino acid as a hotspot, not multiple
  * Open scripts_and_helpers/1_colabfold_sequence_filtering.ipynb
    * Do not run the first cell
    * Run the second cell and wait for the calculation to finish
      * When finished, it will show a histogram that is a visual representation of where the binder wants to bind
        * The reasoning behind this histogram is working with the protein instead of against the protein
        * Using this knowledge, examine the high energy points on the graph. This is where the protein wants to bind
        * If you are seeking a inhibitory protein, find a spot that blocks some sort of active site. This varies widely based on the problem so I cannot give more details
        * If you are seeking a non-inhibitory (more difficult), consider all possible inhibitory actions and find the attractive spot that affects the protein the least. This is a great spot to refer to external resources like youtube/papers to find out what makes a spot inhibitory/non-inhibitory.
    * Using the outputs from the Jupyter notebook, update the hotspots to the site that you determined was best for your needs
    * I would recommend creating a wide range of BINDER_LENGTH so the RFdiffusion has freedom to explore until you learn what is best
    * Run again and see how it goes!

## 5. Locating and Reading the Metrics

All of the outputs are stored in the outputs/ folder in the root directory. When your analysis of a run is complete, I strongly recommend creating an old_results/ folder in the root directory, labeling your run run_1, and moving it into there so you can change the parameters without overwriting your current data. I did not automate this part of the process because I wanted to allow the continued analysis of the current run until you are finished.

To view the metrics, I created a custom script that will analyze and score your creations. Open up scripts_and_helpers/1_colabfold_sequence_filtering.ipynb and run the first cell. As mentioned above, I created an old_results/ folder, but if you don't have that yet, feel free to comment that out, right under the 1. Configuration Section in the first cell (in the SEARCH_DIRS list). 

This will score your protein on the following metrics:

### Metrics Intuition

* ipTM: interface Predicted Template Modeling: This is a score of how confident the interfaces between the chains are in your localcolabfold output. **IMPORTANT NOTE**: If you have more than one chain that is not being designed, it can artificially raise up the ipTM because nature designed a very strong interface between the two or more chains from the wildtype protein.
* pTM: Predicted Template Modeling: This is a score of how confident the model is in the structure stability of the protein. Again, this can also be artificially raised by natural chains, so I often do not use this metric (and it is not used in the scoring cell).
* pLDDT: predicted Local Distance Difference Test: This test is a way to tell how structurally stable the binder is. This is my preferred method of measuring structural stability because it allows for just checking the binder (with some specific code manipulation, done for you in the scoring script, which automatically calculates the binder-specific pLDDT). 
* Rg: A less common metric and one that I implemented myself in the scoring calculation. Rg measures how close the average atom is to the center of mass. Essentially, something like 4 alpha helices close together would have a low Rg (good for my purposes), while a single long straight alpha helix would have a high Rg (bad for my purposes). The reason I implemented this was because I really wanted to get the binder to target a specific location, and I found that the lower Rg genearlly meant both higher specificity (wanting to bind to the hotspots more than anywhere else) and lower clashing (ghosting through the other chains just because it's so long). If Rg is not useful for you, feel free to modify 1_colabfold_sequence_filtering.ipynb to remove it.

### Good Scores for ipTM, pLDDT, and Rg

A frame of reference for the values of each of the scores is provided below:
* ipTM (0-1):
  * Less than 0.5: Pretty terrible
  * 0.5-0.65: Still pretty bad but has the potential to become better
  * 0.65-0.75: Respectable numbers, especially for a de novo protein (created from scratch)
  * 0.75-0.85 Excellent scores, the folding software is convinced that this will bind in real life
  * 0.85-0.9 Incredible scores
  * 0.9-1.0 Suspiciously good, this could either be a natural interface that has incredible binding, or it could be a flaw in your pipeline. If you get above a 0.9, check your work carefully
* pTM:
  * Not dicussed because it is not in the scoring calculation. Feel free to read about it elsewhere and use it if would be of use for you.
* pLDDT (0-100):
  * Less than 60: Weak structure, may not hold together.
  * 60-70: Still fairly weak, but has a bit of potential
  * 70-80: A good spot that means it is likely to hold together in the lab
  * 80-90: Excellent scores, these are very promising
  * 90-95: Absolutely incredible scores, slightly suspicious but possible
  * 95-100: Very suspicious scores, check you work
* Rg (No set range, examples are for protein of around 60 amino acids):
  * Less than 10: Difficult to achieve, very tight protein
  * 10-11: Feasible, very tight protein
  * 11-12: Still a fairly tight protein
  * 12-13: Noticeably less tight, but may still work for you depending on your needs
  * 13+: Very long protein, could have trouble with specificity

## 6. The Iteration Loop

**Read after you complete your first run!**

**Congratulations! You have successfully completed and run and are now analyzing your results!**

Unfortunately, this is slightly difficult to generalize without knowing the specific goals of every person that will use this resource. In general, if you get a score of 60+ and none of the scores are fraudulent and you confirm that it binds where you want it to in PyMOL or other visualizing software, you have successfully formed an excellent binder.

If you get a score of 55+ and it binds where you want it to, that is still an extremely respectable protein.

Under 55, I would recommend running a few more times to see if you can get it to increase.

As a note, be warned that building de novo proteins can be extremely difficult. Learning how to use these tools together and build this pipeline took me hundreds of hours. Competing in the competition took me thousands of hours of compute time on an RTX 5070 OC to get 1 score over 60, and 3 scores between 55 and 60. Some of that time was wasted while I was learning how to do things properly, but it also just takes a huge amount of iterations to find viable binders, which is why much of drug development is an open problem, despite incredible recent advances in the field.

Here are few notes for your next steps:
* If you have a low ipTM but high pLDDT (or the protein is binding in the wrong spot), the binder likes its structure, but does not like the interface to the protein. Try running the histogram again and finding a new hotspot. You could also try changing the specific hotspots that you listed in the config. For example, instead of "A1,A2,A3,A4,A5,A6" when you see a spike near A3, try something like "A1,A4,A6", where those specific hotspots are available for binding. I would recommend using sticks mode in PyMOL to help you determine which should be used at hotspots.
  * In PyMOL, use *show sticks* on the target protein's surface. Look for pockets or clusters of amino acids that look like they could grip a binder. Those are your best candidates for HOTSPOTS.
* If you have a high ipTM but low pLDDT, the binder is very happy with the interface, but does not believe it will be structurally stable. One way to fix this is to increase the MPNN_TEMP to 0.3 and the SEQ_PER_BACKBONE to 20. This will generate more fasta sequences and choose the more structurally sound of them all, though it might slightly worsen the ipTM as the code prioritizes structural stability.
* If either or both scores are bad, you could also consider modifying the binder length. If it is too low or high, both ipTM and pLDDT might suffer. This can be found through experimentation.

**You will see 99% of your binders fail. This is expected. Just modify a few parameters and run it again!**


## 7. The PyRosetta Toolkit

**Congratulations on Finding a Successful Binder!** Now, it is time to optimize. PyRosetta is amazing for protein mutation and analysis. Open the 2_pyrosetta_general_analysis.ipynb.

While the cells are labeled and commented, I will give extra clarification on any steps that I feel may be confusing.

Cell 3: Calculates the total system energy of your protein and binder. A lower score is better. The goal is to get into the far negatives, but where exactly depends on the problem. You will likely see a very high number in cell 3. This is perfectly okay as the folding software put sidechains just a bit too close together, and now they are all repelling. Luckily, a quick FastRelax takes care of this.

Cell 4: FastRelax: This will gently wiggle the structure until the energy is lowered. This should drastically reduce the values from Cell 3. Another incredibly important metric to utilize: ddG (delta delta G). This means the change in Gibbs Free Energy when the binder is connected to the protein vs when they are forcibly separated. A more negative value is better. Generally, the goal is around -30 REU (Rosetta Energy Units), but it can vary. Intuitively, this represents how much the binder wants to stick to the protein.

**NOTE**: If your ddG is in the millions, you likely have a disulfide bond. This depends on your problem, but is has a potential to be bad, so do some research and consider mutating it away.

Cell 5: This will calculate the binding energy and find the top 5 most stabilizing and destabilizing residues (amino acids). 

Cell 6: An interactive viewer will allow you to move the protein around and look at how it bonds.

Cell 7: This is when the power of Rosetta comes to the fore. PackRotamersMover is an incredibly powerful function that automatically mutates the protein, relaxes it, and checks if it is better. It repeatedly does this until it reaches a local energy minima, very often finding a better setup for your binder. This is the slowest part of PyRosetta, but it is very worth it for all of your high-quality binders from 1_colabfold_sequence_filtering.ipynb. Later, I might set up a pipeline to automate this part for scores above a certain threshold, but I have not done this yet.

Cell 9: Essentially just cleans and outputs the results from PackRotamersMover. As of now, this is still hardcoded to nitrogenase, but I may try to fix that in the future. Just be aware that other chains in the output are likely incorrect for your protein, but the final chain is accurate.

## 8. Utilizing AlphaFold3

**Most Useful After PyRosetta**: Now, you have high quality binders that have been further optimized by PyRosetta. Unfortunately, due to the changes in the fasta sequence, the structure previously predicted by LocalColabFold may no longer be accurate, so you need to re-run a folding software. This is when AlphaFold3 is very useful, particularly if you needed to crop the protein. Extract the full fasta sequence from your specific protein, add on another chain with your custom fasta that was the output from PyRosetta and run it through AlphaFold3 server (free for 30 proteins per day, which is a good number). Because this step is only intended for the final candidates, you should not need to run AlphaFold3 more than 50 times. As such, it has not yet been automated. This may change in the future. If you did not need to crop your protein and are happy with running it all locally, feel free to copy, paste, and modify 4_run_localcolabfold.sh to run on your folder of pyrosetta_outputs. AlphaFold3 has the ability to handle thousands of amino acids in just a few hours, which makes it an incredible tool for final processing. It also may have higher fidelity. 

Once you have a final result from AlphaFold3 with high scores and binds in the correct location, download it to your computer. I would recommend opening all 5 cifs in PyMOL (PyMOL can open cifs just like it can open pdbs) and aligning them to check if all of the models show your bindering attaching in the same spot. If so, success! You have developed an awesome new protein. If not, it may still be okay with a high magnitude (very negative) ddG on the original spot, or if the other locations are also inhibitory/non-inhibitory, based on your needs. Well done!

## 9. Conclusions

**Thank you for taking the time to utilize this repository!** I truly hope you learned a lot and were able to create something amazing. Please feel free to reach out to discuss this repository at tony_peonio@mines.edu.

## Acknowledgements

**Technical Foundations**
This pipeline is a wrapper that integrates several groundbreaking tools from the computational biology community. This work would not be possible without
* **[RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)**: For the generative design of protein backbones.
* **[ProteinMPNN](https://github.com/dauparas/ProteinMPNN)**: For fast and accurate sequence design.
* **[LocalColabFold](https://github.com/YoshitakaMo/localcolabfold)**: For local, high-throughput structural prediction and validation using AlphaFold2.
* **[AlphaFold3](https://isomorphiclabs.com/alphafold3)**: For final high-fidelity validation.
* **[PyRosetta](https://www.pyrosetta.org/)**: For chemical interface analysis and energy scoring.

**Personal Acknowledgements**
A special thank you to those who helped shape this project from a competition pipeline into a public resource:
* Abby Miller: For her help as a Quantitative Biological Engineer, ensuring the facts were straight and the explanations were clear to a scientific audience. She also helped by teaching me the biology behind these tools as I created this software.
* Jessica Peonio: For providing the "layman's perspective", catching typos and increasing clarity for those outside the field of biology.