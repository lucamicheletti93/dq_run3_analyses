## pp reference cross section

- Script to compute the cross section:
```ruby
  python computeXsec.py configs/config_xsec_charmonia.yml --run
  ```
- Script to compute the Acceptance x Efficiency:
```ruby
  python computeAxe.py configs/config_axe_charmonia.yml --run
  ```
- Script to perform the interpolation procedure:
```ruby
  python interpolation.py configs/config_interpolation.yml --run
  ```
