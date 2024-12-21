# C++ Code Minimum Rooted Drone Deployment Problem with Neighborhoods (MDPN)

![C++](https://img.shields.io/badge/-C++-blue?logo=cplusplus)
![License](https://img.shields.io/badge/license-MIT-green)

## Overview

This repository contains the C++ implementation of the algorithms and experiments described in our paper, **"Single- and Multi-Depot Optimization for UAV-Based IoT Data Collection in Neighborhoods"**, which has been accepted for publication in **ACM Transactions on Sensor Networks**.

In this paper, we investigate the problem of deploying the minimum number of Unmanned Aerial Vehicles (UAVs) and determining their flying tours to collect data from all Internet of Things (IoT) sensors. We study this problem in a neighborhood scenario where a UAV can collect data from an IoT sensor if the distance between them is less than the wireless communication range of the IoT sensor. Since UAVs are powered by batteries with limited energy, we assume that the total energy consumed during the flying tour of each UAV is bounded by a given budget.

We present the **Minimum Rooted Drone Deployment Problem with Neighborhoods (MDPN)**, which is NP-hard, and propose two approximation algorithms for the single-depot case. One of these algorithms is a bi-criteria approximation that guarantees a solution where the tour's cost is violated by a factor of \(1 + \epsilon\). Furthermore, we extend these two algorithms to the multi-depot scenario. Finally, we evaluate our algorithms in three different scenarios: an ideal one where the communication range is a circle and the data transfer rate is constant, and two more realistic scenarios where we introduce irregularity in the communication range and a non-constant data transfer rate.

---

## Installation

Clone the repository and follow these steps to build the project:

```bash
git clone https://github.com/bettisfr/MDPN_TOSN.git
cd MDPN_TOSN
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/TOSN
```

---

## Citation

If you use this repository, please cite our paper:

```bibtex
@article{betti2025TOSN,
  title={Single- and Multi-Depot Optimization for UAV-Based IoT Data Collection in Neighborhoods},
  author={Betti Sorbelli, Francesco and Ghobadi, Sajjad and Pinotti, Cristina M.},
  journal={ACM Transactions on Sensor Networks},
  year={2025},
  issn = {1550-4859},
  url = {https://doi.org/10.1145/3704810},
  doi = {10.1145/3704810},
  note={To appear}
}
```

---

## Contact Us

For questions, feedback, or collaboration opportunities, feel free to reach out to us:

- **Francesco Betti Sorbelli**: [francesco.bettisorbelli@unipg.it](mailto:francesco.bettisorbelli@unipg.it)
- **Sajjad Ghobadi**: [sajjad.ghobadibabi@unipg.it](mailto:sajjad.ghobadibabi@unipg.it)
- **Cristina M. Pinotti**: [cristina.pinotti@unipg.it](mailto:cristina.pinotti@unipg.it)
