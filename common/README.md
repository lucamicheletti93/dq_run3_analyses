## Common tools
### Tracking efficiency
- Create the histograms from the output of the mch tracking efficiency task:
  ```ruby
  python trackingEfficiency.py configs/config_tracking_efficiency.yml --hist
  ```
- Compute the tracking efficiency in data and MC and compare:
  ```ruby
  python trackingEfficiency.py configs/config_tracking_efficiency.yml --run
  ```
