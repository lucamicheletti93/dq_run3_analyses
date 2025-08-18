# Simulation code: how to use

Run the simulation locally:

- Load O2Sim
- Run the code:
  ```ruby
  ./runOniaSim.sh
  ```

- WARNING: check always the number of events you want to produce before running
- WARNING: do not push the produced file, it will overwrite the existing one!!!

Process the output of the simulation:

```ruby
root -l readSim.c
```