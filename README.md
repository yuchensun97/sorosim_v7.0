# Sorosim 7.0

This is a `MATLAB` implementation of the **Real-time dynamic model of Soft Manipulators with Cross-section Inflation**. The source code is modified from the original work of [SoRoSim](https://github.com/Ikhlas-Ben-Hmida/SoRoSim) to include the cross-sectional inflation effect of the classic Cosserat rod model using the Geometric Variable Strain (GVS) approach.

We provide three examples to demonstrate the capability of the model, which can be found in the `./Custon` folder. The examples include the stiffness tuning of a soft manipulator, the dynamic simulation of a octopus'arm reaching movement, and the fetching motion.

The detailed description of the model are available from the author upon reasonable request.

If you wish to use this code for your research, please properly cite the following paper.

```
@article{sun2025real,
  title={Real-time dynamics of soft manipulators with cross-sectional inflation: application to the octopus muscular hydrostat},
  author={Sun, Yuchen and Mathew, Anup Teejo and Afgan, Imran and Renda, Federico and Laschi, Cecilia},
  journal={Proceedings of the Royal Society A},
  volume={481},
  number={2314},
  pages={20240642},
  year={2025},
  publisher={The Royal Society}
}

@ARTICLE{9895355,
  author={Mathew, Anup Teejo and Hmida, Ikhlas Ben and Armanini, Costanza and Boyer, Frederic and Renda, Federico},
  journal={IEEE Robotics & Automation Magazine}, 
  title={SoRoSim: A MATLAB Toolbox for Hybrid Rigid–Soft Robots Based on the Geometric Variable-Strain Approach}, 
  year={2023},
  volume={30},
  number={3},
  pages={106-122},
  keywords={Robots;Mathematical models;Strain;Soft robotics;Robot kinematics;Matlab;Analytical models},
  doi={10.1109/MRA.2022.3202488}}
```

## License
This project is licensed under the [MIT License](LICENSE).
