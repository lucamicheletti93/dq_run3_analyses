## Common tools
### Tracking efficiency
- Run the corresponding task locally after modifying and copiling the task tableMakerMuonMchTrkEfficiency (temporary fix: comment all reserve command in the task):
  ```ruby
  ./run_train.sh --skip-perf
  ```
- Create the histograms from the output of the mch tracking efficiency task:
  ```ruby
  python trackingEfficiency.py configs/config_tracking_efficiency.yml --hist
  ```
- Compute the tracking efficiency in data and MC and compare:
  ```ruby
  python trackingEfficiency.py configs/config_tracking_efficiency.yml --run
  ```
