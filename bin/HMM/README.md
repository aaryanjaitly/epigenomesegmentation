# C++ Implementation of a flexible distribution HMM for chromatin segmentation

Implementation of a multivariate HMM, supporting several distributions (**PO, BI, NBI, BB, BNB, SI, ZAP, ZANBI, ZASI, ZABNB**) to calculate the emission probability in a state for discrete data, like histone modifications and optionally DNA methylation as additional input (two columns, distributions: **BI, BB, AB**). 

## Dependencies

- C++ compiler with support for C++ 17. 
- CMake version >= 3.16 (<https://cmake.org/>)
- Boost >= 1.55.0 (<https://www.boost.org/>)
- OpenMP (on Linux if not available install with `sudo apt-get install libomp-dev`)
- ROOT (<https://root.cern/>)

## Compilation

```
mkdir build
cd build
cmake ..
make -j 4
```

To also build the tests use instead 
```
cmake .. -DGTEST_SRC_DIR=/usr/src/googletest/
```
depending on the path to googletest. The tests can then be run with `make test`.

## Parameters
### Executable **./HMMChromSeg**
```
  -h [ --help ]                         Output help message.
  -m [ --modelInput ] arg               Serialized HMM with starting 
                                        parameters.
  -r [ --regions ] arg                  Start indices for independent regions 
                                        (starting at zero and whitespace 
                                        separated), if not set no splitting is 
                                        performed.
  -c [ --countMatrix ] arg              Count matrix for training or decoding.
  -x [ --methylationMatrix ] arg        Optional matrix for DNA methylation, 
                                        only required if model has set 
                                        methylation to true.
  -t [ --train ]                        Should the HMM model be trained first?
  -i [ --maxIteration ] arg (=300)      Maximum number of iteration during 
                                        training.
  -j [ --iterationsPretrain ] arg (=5)  Number of iterations before topology is
                                        adjusted.
  -n [ --numAdjustments ] arg (=2)      How often should the topology be 
                                        adjusted?
  -e [ --epsilon ] arg (=0.1)           Convergence of likelihoods as  
                                        termination criterion during training.
  -o [ --modelOutput ] arg              Optional output file for parameters of 
                                        new HMM (only required if HMM is 
                                        trained first).
  -v [ --viterbiPath ] arg              Optional output file for state sequence
                                        using Viterbi decoding.
  -d [ --posteriorDecoding ] arg        Optional output file for state sequence
                                        using posterior decoding.
  -p [ --threads ] arg (=1)             Number of threads for training the HMM 
                                        model.
```

### Executable **./TopologyHMM**

```
Command line options:
  -h [ --help ]                         Output help message.
  -m [ --modelInput ] arg               Serialized HMM with starting 
                                        parameters.
  -r [ --regions ] arg                  Start indices for independent regions 
                                        (starting at zero and whitespace 
                                        separated), if not set no splitting is 
                                        performed.
  -c [ --countMatrix ] arg              Count matrix for training or decoding.
  -x [ --methylationMatrix ] arg        Optional matrix for DNA methylation, 
                                        only required if model has set 
                                        methylation to true.
  -t [ --train ]                        Should the HMM model be trained first?
  -a [ --adjustTopology ]               Should the HMM topology be updated 
                                        during training?
  -i [ --maxIteration ] arg (=300)      Maximum number of iteration during 
                                        training.
  -j [ --iterationsPretrain ] arg (=5)  Number of iterations before topolgy is 
                                        adjusted.
  -e [ --epsilon ] arg (=0.10000000000000001)
                                        Convergence of likelihoods as 
                                        termination criterion during training.
  -o [ --modelOutput ] arg              Optional output file for parameters of 
                                        new HMM (only required if HMM is 
                                        trained first).
  -v [ --viterbiPath ] arg              Optional output file for state sequence
                                        using Viterbi decoding.
  -d [ --posteriorDecoding ] arg        Optional output file for state sequence
                                        using posterior decoding.
  -p [ --threads ] arg (=1)             Number of threads for training the HMM 
                                        model.
```

## Example Input

- count matrix (tab separated)
```
0   0   0
12  26  30
14  24  28
...
```

- count matrix methylation (tab separated, first column corresponds to coverage and second the methylation)
```
5	4
8	7
1	1
0	0
0	0
...
```

- model input (json file, for classical model)

```json
{
    "states": "4",
    "marker": [
        "H3K36me3",
        "H3K9me3",
        "H3K27me3"
    ],
    "emission": [
        [
            {
                "distribution": "PO",
                "parameters": {
                    "lambda": "82.882550792173276"
                }
            },
            {
                "distribution": "PO",
                "parameters": {
                    "lambda": "12.219109112889196"
                }
            },
            {
                "distribution": "PO",
                "parameters": {
                    "lambda": "8.4965508552604874"
                }
            }
        ],
        ...
    ],
    "initial_state_distribution": [
        "0.25",
        "0.25",
        "0.25",
        "0.25"
    ],
    "transition_matrix": [
        [
            "0.87410650925151345",
            "0.12351866650378064",
            "0.00020746336180438969",
            "0.0021673609348816093"
        ],
        ...
    ]
}
```
If **initial_state_distribution** or **transition_matrix** are not specified, they are initialized uniformly.

- model input (json file, for topology model)

```json
{
    "states": "4",
    "topology": {
        "1": [
            "0",
            "1"
        ],
        "2": [
            "2"
        ],
        "3": [
            "3",
            "4",
            "5"
        ],
        "4": [
            "6"
        ]
    },
    "marker": [
        "H3K36me3",
        "H3K9me3",
        "H3K27me3"
    ],
    "emission": [
        [
            {
                "distribution": "NBI",
                "parameters": {
                    "p": "0.15300346425732353",
                    "r": "0.45033512558679911"
                }
            },
            {
                "distribution": "BNB",
                "parameters": {
                    "alpha": "61610.675832895198",
                    "beta": "188636.48043647193",
                    "r": "2.7396825337719624"
                }
            },
            {
                "distribution": "SI",
                "parameters": {
                    "mu": "0.0022781493817052301",
                    "sigma": "231.58815040847352",
                    "v": "1.3480227804831695"
                }
            }
        ], ...
    ],
    "initial_state_distribution": [
        "0.25",
        "0",
        "0.25",
        "0.25",
        "0",
        "0",
        "0.25",
    ],
    "transition_matrix": [
        [
            "0.9263",
            "0.0736",
            "0",
            "0",
            "0",
            "0",
            "0"
        ],
        [
            "0",
            "0.9263",
            "0.0574",
            "0.0093",
            "0",
            "0",
            "0.007"
        ],
        ...
}
```