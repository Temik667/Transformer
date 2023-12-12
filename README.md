# Transformer
Here you can see the code for the transformer winding monitoring using sweep frequency response analysis and machine learning.
The key reason for using machine learning in this scenario is that people often rely on the expert's knowledge when analysing the damages in the transformers.
Therefore this imposes an ambiguity and sometimes might lead to irreversible consequences costing significant financial losses. 
Therefore, our project aims to eliminate the human factor and propose a reliable analysis method for the industry.

# Structure
1. **Sample generators** are the models written in Python to be used to simulate different damages. The values are swept to create data for training an ML model.
2. **Netlists** are LTSpice code(simulation software for electrical circuits) that represents transformer models. These are then run to get data
3. **Netlists_data** is a folder with the frequency responses of all possible damages. A plotted example can be seen in the Healthy_transformer.jpg
4. **ML** is the folder with all the training
