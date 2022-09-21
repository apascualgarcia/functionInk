---
layout: page
title: Visualization
subtitle: Last update September 2022
---

# Visualization

In this page we explain how to represent the communities found with functionink. We will use the results obtained in the [Vignette](../Vignette) which can also be found in the directory `vignette_example`.

There are many methods to plot networks. Here we illustrate a procedure with the popular software [Cytoscape](https://cytoscape.org/) (version 3.3.0).

First, we start Cytoscape and we select the option "From network file":

![ Figure 1](../../_images/visualization/Cytoscape1.png)

We load the file `Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.format2.txt` located in the directory `data`. In "Advanced Options" we add "ignore lines starting with #" to skip the header, and unclick the option "Use first line as column headers" (_Note: Cytoscape sometimes does not react to these options, just click and unclick a couple of times and you will have it_)

![ Figure 2](../../_images/visualization/Cytoscape2.png)

Then, we edit the fields Column 1 to Column 4, indicating which is the "source node" (Column 1), the "target node" (Column 2), and fixing the type of entry to "edge attributes" for columns 3 and 4 giving a name to Column 3 (interaction) and Column 4 (type):

![ Figure 3](../../_images/visualization/Cytoscape3.png)

After renaming the columns, it should look something like this:

![ Figure 4](../../_images/visualization/Cytoscape4.png)

The network is loaded, and it looks pretty messy. Now we load the partition we obtained with functionink to shed some light on the structure. We will load the file `Partition-NL_Average_StopStep-47_Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt` found in the directory `vignette_example`. To load it we do File > Import > Table > File and a prompt window opens. Again, we go to Advanced Options and we proceed as we did with the network (include the hash to ignore and include the first line). After editing the columns' names it should look like this:

![ Figure 5](../../_images/visualization/Cytoscape5.png)

We are now ready to identify the communities. To simplify the visualization, we can remove competitive interactions within each pool since were randomly drawn and hence we may want to explore first mutualistic interactions, which have a high nestedness. We filter those edges with "type" equal to 2, which identify competitive links. In the tab "Select" of the Control Panel, you  can just  add one condition to filter them: 

![ Figure 6](../../_images/visualization/Cytoscape6.png)

After "apply" the filter, in View > Edges > Hide selected edges, you will hide them. And we finally identify the communities. We go to layout > Group attributes layout > functionink, and Cytoscape will separate the communities.

![ Figure 7](../../_images/visualization/Cytoscape7.png)

From this representation, with a bit of manual reordering and the aesthetic possibilities that Cytoscape offers you can make a cleaner figure, as the one shown in the publication of the method. For example, in this figure we coloured the nodes according to their community, so if they are close in space and have the same colour they belong to the same community

![ Figure 8](../../_images/visualization/Cytoscape8.png)

where it is apparent the definition of guilds, except perhaps for those nodes with only one link that sometimes belong to the same community and sometimes don't. The reason is that their split depends on the competitive interactions. So let's bring them back and colour them differently.

![ Figure 9](../../_images/visualization/Cytoscape9.png)

And indeed these communities are split because they have different competitive links. You may want to explore more in detail specific nodes to verify the definition of guild. The Cytoscape project is found in the folder `vignette_example/figures` with the name `Vignette_Cytoscape.cys`

