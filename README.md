# gps-lab3

gps-lab3 is a repository of MATLAB code associated with Lab 3 of Fundamentals of GPS.

## Installation

Clone the repo.

```sh
git clone git@github.com:tannerkoza/gps-lab3.git
```

## Usage

Once cloned, open the project as the current folder in MATLAB and run the following in the Command Window:

```sh
startup
```
This adds all subfolders in the project to the MATLAB path. This same functionality is also automatically achieved when MATLAB starts in this project folder.

If other receiver data is needed, add the unparsed receiver `.mat` files to `data/receiver_data`. Then, run `lab_parser` in the Command Window from the root folder. The parsed files will be saved in `data/parsed_data` to be used in other scripts.  
## Contributing
Contributions are PROHIBITED as this is a graded assignment. 
## License
This repo serves as a future reference to concepts learned in this class and as personal repo management practice. It is not meant to be cloned and turned in.

[MIT](https://choosealicense.com/licenses/mit/)
