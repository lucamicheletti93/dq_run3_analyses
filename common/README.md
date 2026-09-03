## Common tools
### Tracking efficiency
- temporary fix: comment all `reserve` command in the task `tableMakerMuonMchTrkEfficiency`
- Run following script in the reposotories `/task_config_mch_trk_eff_data` and `/task_config_mch_trk_eff_mc`:
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
