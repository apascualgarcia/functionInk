---
layout: page
title: Issues
subtitle: Last update September 2022
---

## Frequent problems

* Your file is not formatted correctly:
   * It is not tab-separated.
   * The header does not start with `#`.
   * The columns do not have the correct order.

* Your file has weights equal to zero. 
   * The links does not exist, and the algorithm does not know how to deal with it, these links should be removed.

* The algorith returns a message like `The value for network weight argument (-w) is not numeric = 1`. In that case either:
   * You do not have installed the module `Scalar::Util`
   * You are working in an environment with a Perl version that do not have the necessary modules installed (see page [Install](../Install)).

## Other issues

If your problem is not listed, you may want  to open an [issue](https://github.com/apascualgarcia/functionInk/issues) and I will try to help you!
