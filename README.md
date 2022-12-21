# propsim
A simple and intuitive tool for simulating different types of aircraft engines.

##  Overview
This library aims to bring a simple and intuitive way to perform computational calculations for the design and validation of the main models of aeronautical engines.

Through a core `AircraftEngines` class the user is able to switch between the following engine types and evaluation methodologies:

Engine Model  | Functionality                                              | Implemented function
------------- | -------------                                              | -------------
Turbojet      | Ideal<br />Non-ideal on design<br />Non-ideal off design   | `ideal_turbojet`<br />`real_turbojet`<br />`real_turbojet_off_design` 
Turbofan      | Ideal<br />Non-ideal on design<br />Non-ideal off design   | `ideal_turbofan`<br />`real_turbofan`<br />`real_turbofan_off_design` 
Turboprop     | Non-ideal on design<br />Non-ideal off design              | `ideal_turboprop`<br />`real_turboprop`
Ramjet        | Ideal                                                      | `ideal_ramjet`

##  Setup process

First, one must install the library using the following command in an environment containing Python 3.6 or higher.

```
pip install propsim
```
This command will install the library and its dependencies.

##  Usage cases

Once the library is properly installed, we have the following use cases:

1. The first one consists of carrying out a single point analysis with a fixed $\pi_c$ value. In this specific example an ideal turbofan type engine is shown.

    ```
    from propsim import AircraftEngines

    engines = AircraftEngines(12500)

    engines.ideal_turbofan(M0=0.7, gamma=1.4, cp=1004, hpr=42.8*10**6, Tt4=1850, pi_c=10, pi_f=2, alpha=5)
    ```

    The expected result for this test case is the following data set:

    ```
    {
    'pi_c': [10],
    'F_m0': [279.62], 
    'f': [0.03], 
    'S': [1.94e-05], 
    'eta_T': [0.52], 
    'eta_P': [0.46], 
    'eta_Total': [0.24], 
    'FR': [4.46]
    }
    ```


2. The second one consists of carrying out a bach analysis were $\pi_c$ values vary and so the output data can be plotted for carrying out the appropriate analysis. This feature is currently available only for the ideal turbofan, ideal turbojet and real turbojet. In this specific example an ideal turbofan type engine is also shown.

    ```
    from propsim import AircraftEngines

    engines = AircraftEngines(12500)

    engines.ideal_turbofan(M0=0.7, gamma=1.4, cp=1004, hpr=42.8*10**6, Tt4=1850, pi_c=10, pi_f=2, alpha=5, batch_size=3, min_pi_c=7, max_pi_c=15)
    ```

    The expected result for this test case is the following data set:

    ```
    {
    'pi_c': [7, 9.66, 12.33, 14.99], 
    'F_m0': [271.16, 278.94, 283.34, 286.05], 
    'f': [0.031, 0.032, 0.0312, 0.031],
    'S': [2.06e-05, 1.95e-05, 1.87e-05, 1.82e-05],
    'eta_T': [0.47, 0.52, 0.55, 0.57],
    'eta_P': [0.48, 0.47, 0.46, 0.45],
    'eta_Total': [0.23, 0.24, 0.25, 0.26],
    'FR': [4.17, 4.43, 4.58, 4.67]
    }
    ```