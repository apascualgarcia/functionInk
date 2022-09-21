---
layout: page
title: Install
subtitle: Last update September 2022
---

### Requirements

The method uses two scripts written in R and one script  in R. If you haven't used Perl in the past no worries, it is a language similar to R and you just need to have a Perl interpreter installed in your computer. Running the scripts is then easy and do not require any further installation. Similarly, the script for the analysis requires R.

Most Unix distributions come with a Perl interpreter so you will possibly not need to install anything, if you are running a different OS you can find more information [in this page](https://perldoc.perl.org/5.32.0/perlfaq2.html#What-machines-support-Perl%3f-Where-do-I-get-it%3f) on how to install Perl in your comp.

The script uses two modules (`POSIX` and `Scalar::Util`) and the latter may not come installed alongside your Perl native installation. An easy way to install modules is via [CPAN](https://www.cpan.org/modules/INSTALL.html). Once you have CPAN installed you can install modules (see [this post](https://perlmaven.com/how-to-install-a-perl-module-from-cpan) for more details and installation in operating systems different than Linux and [this post](https://ostechnix.com/how-to-install-perl-modules-on-linux/) for other alternatives in Linux).

We show how to install the module required. In Linux, first open a CPAN shell:

```
$> sudo perl -MCPAN -e shell
```

And then install, the required module:

```
$CPAN> install Scalar::Util
```

### Clone the repo

Once you installed the modules described above and the Perl interpreter, you could simply download the scripts mentioned in the main page [overview](/) from the [repository](https://github.com/apascualgarcia/functionInk). If you want to follow the [Vignette](../Vignette), it is easier to directly clone the repository in your computer:

```
$> git clone https://github.com/apascualgarcia/functionInk.git
```

We then go to the root of the repository:

```
$> cd path_to_the_repository
```

And we give permission to execute the scripts:

```
$> chmod 755 NodeLinkage.pl NodeSimilarity.pl 
```

### Test the scripts

We start testing the script to create the similarity matrix between nodes. In the directory `fake_example` we
have a small network for testing. We run (from the root of the directory):

```
./NodeSimilarity -w 1 -d 1 -t 1 -f fake_example/network1.tsv
```

After running, you should find the file ```Nodes-Similarities_network1.tsv``` in your folder, containing the similarities. We now test the agglomerative clustering:

```
./NodeLinkage.pl -fs Nodes-Similarities_network1.tsv -fn fake_example/network1.tsv
```

with no options it creates the files ```HistExtend-NL_Average_NoStop_network1.tsv``` and ```HistCompact-NL_Average_NoStop_network1.tsv```. Now we test that it works with more options, we make it stop at step 2 of the clustering:

```
./NodeLinkage.pl -fs Nodes-Similarities_network1.tsv -fn fake_example/network1.tsv -s step -v 2
```

And you will find four more files with the label ```StopStep-2```. Now you may want to go to the [Vignette](../Vignette) for a more spicy example with further explanations.



