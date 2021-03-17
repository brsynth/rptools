# rpThermo

Calculate the formation energy of chemical species and the gibbs free energy of reactions and the heterologous pathway. This tool uses the [component contribution](https://gitlab.com/elad.noor/component-contribution) method for determining the formation energy of chemical species that are either not annotated, or cannot be found in the internal database.  

## Input

Required:
* **-input**: (string) Path to the input file
* **-input_format**: (string) Valid options: tar, sbml. Format of the input file

Advanced Options:
* **-pathway_id**: (string, default=rp_pathway) ID of the heterologous pathway

## Output

* **-output**: (string) Path to the output file 

## Dependencies

* [Marvin:](https://chemaxon.com/products/marvin)

## Building the docker

NOTE: you need to have a valid [Marvin](https://chemaxon.com/products/marvin/download) account and Marvin licence (named license.cxl) in the root directory. Furthermore, the Dockerfile needs to be modified to have the addition of the source.list as per the deb instructions in the official Marvin website.

```
docker build -t brsynth/rpthermo-standalone -f Dockerfile .
```

## Contributing

TODO

## Versioning

v0.1

## Authors

* **Melchior du Lac** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

### How to cite rpThermodynamics?

TODO