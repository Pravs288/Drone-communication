# Supplementary Material to: Location Optimization for Tethered Aerial Base Station Serving mmWave High Altitude UAVs.

### Pravallika Katragunta, Carleton University, School of Computer Science, Canada.

### Michel Barbeau, Carleton University, School of Computer Science, Canada.

### Joaquin Garcia-Alfaro, Samovar, Institut Polytechnique de Paris, Telecom SudParis, 91120 Palaiseau, France

### Evangelos Kranakis, Carleton University, School of Computer Science, Canada.

### Venkata Srinivas Kothapalli, Ericsson Canada, K2K 2V6, Ottawa, Ontario, Canada.

## Summary
Uncrewed Aerial Vehicle-User Equipment (UAV-UE) is integral to millimeter wave (mmWave)-based wireless cellular systems. 
UAV-UE at high altitudes encounter limited connectiv- ity with terrestrial base stations. Tethered Aerial Base 
Stations (TABS) are viable alternatives to terrestrial base stations. Optimal placement of a TABS in a three-dimensional 
environment is necessary and critical to serve multiple moving UAV-UE units with reliable connectivity. In this work, 
we propose a contextual multi- armed bandit framework to learn the optimal TABS locations. We consider multiple UAV-UE 
units moving at high altitudes in an uplink mmWave setting. Under this framework, the TABS acts as a learning agent 
leveraging position information about served UAV-UE units to provide connectivity with minimum Signal to Noise 
Ratio (SNR) threshold requirements. We first compare the Upper Confidence Bound (UCB) and Thompson Sampling (TS)-based 
learning strategies against the traditional naive-based approach. Our simulation results show that the TS-based approach 
learns optimal locations with a 31% and 51% average regret-reduction ratio (ARR) over UCB and naive-based approaches, 
respectively. Also, the TS-based learning strategy for TABS reliably achieves the required SNR for UAV-UE units under 
multiple contexts, compared to a static TABS location.


## Reference

If using this code for research purposes, please cite:

P. Katragunta, M. Barbeau, J. Garcia-Alfaro, E. Kranakis, V. S. Kothapalli. 
Location Optimization for Tethered Aerial Base Station Serving mmWave High 
Altitude UAVs. 2024 IEEE Canadian Conference on Electrical and Computer 
Engineering (CCECE 2024), Kingston, Ontario, Canada, August 06–09, 2024.

```
@inproceedings{katragunta2024location,
  title={{Location Optimization for Tethered Aerial Base Station Serving mmWave High Altitude UAVs}},
  author={Katragunta, Pravallika and Barbeau, Michel and Garcia-Alfaro, Joaquin and Kranakis, Evangelos and Kothapalli, Venkata Srinivas},
  booktitle={2024 IEEE Canadian Conference on Electrical and Computer Engineering (CCECE 2024)},
  pages={271--276},
  year={2024},
  organization={IEEE},
  doi={10.1109/CCECE59415.2024.10667117},
  url={https://doi.org/10.1109/CCECE59415.2024.10667117},
}
```


