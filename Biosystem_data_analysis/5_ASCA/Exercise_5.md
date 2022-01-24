###Q4: Compare the variance of the individual variables with and without taking the square root of the data.

### Q3: Why is it important for PCA that the measured values of all variables are approximately of the same magnitude?

### Q1:

### · How many glucosinolates are measured?

11 glucosinolates

### · What are the levels of the factor *Time?*

4

### · What are the levels of the factor *Method?*

3

### · How many samples were taken for each *Time/Method* combination?

60

### Q2: Discuss the problematic issues, specific for metabolomics data that you observe for this data.

*We see large differences in intensity between the variables. There are also skewed distributions for variables GBC, NEO. Thus there are few samples with very large intensity values, while most of the measurements are below the average value. This needs to be corrected for before a PCA is performed.* **(From Answer)**

### From the box plot you see that some variables clearly show the effect of the treatment while others seem unaffected. Which variables show a clear treatment effect?

PRO,  RAPH, GNL, X4OH, GBN, X4MeOH, 



### Q3: Why is it important for PCA that the measured values of all variables are approximately of the same magnitude?

*As  PCA  focusses  to  describe  variation  in  the  data,  it  prefers  columns  with  much 
variation.  If  all  variables  have  equal  variance  (which  is  the  case  after  standardization 
(also called autoscaling)) then all variables have an equal chance of being included in 
the model. In that case the correlation of the variable with all other variables decides 
its contribution to the model.* **(From Answer)**

### Why is centering a good idea before doing PCA?

*As we are interested in variation between samples and not how samples deviate from 0. 
The centering step places the origin in the center of the data, so PCA concentration on the 
direction of optimal variation between the samples.* **(From Answer)**

### Do PCA on the data without centering. How can you see from the score plot that the data are not centered?

*0 is not in the center of the data on the first PC, because data is not centered. Note  
the 2nd PC is centered. This comes from the orthogonal restriction of the PC’s.* **(From Answer)**

### Which variables contribute most to the variation? If this surprises you, have another look at the box plots.

From the loadings plot, NEO, GBC

From the box plot, we can observe that NEO and GBC has the largest difference on mean between treatment groups. So this might be the reason.



### Q4: Compare the variance of the individual variables with and without taking the square root of the data.

### Also compare the PCA score and loading plots (biplot) for data with and without taking the square root for both centered and non-centered data.

### What is your conclusion?

By taking the square root of the data, we reduced the difference of variance between different variable groups, although GBC and NEO still have the largest variance.

After taking square root, the variation contributed by variables other than GBC and NEO has increased.



*Transformation removed the skewed distributions, Centering puts the center in the middle of the data to better find the optimal directions between the samples.*

***Transformations also have some scaling effect making the variation for variables more similar.***

### Q5:

skip **(READ THE ANSWER!)**

NOTE: ***In ASCA, the PCA analysis is performed on group averages!!!!!!***

### Q6: Which factor explains most of the variation?

factor 2, the treatment group

### How much of the variation is in the residuals?

40.35%

### Check how the values change when the data are not centered after square root transformation.

The overall means contribute to most variation, while factor 1, 2 and residuals become trivial.



### Q7: How much variation is explained by each of the factors?

10.22 for factor 1, 49.42 for factor 2

### Seeing the score plots, what would you say about the effect of the two treatments (time and application)?

The effect of time is relatively small comparing with the effect of treatment. In the score plot, the clustering is very obvious for treatment comparing with time group.

### Compare the results you find with ASCA with those you found by observing the box plots. Especially pay attention to the treatment effect on the different glucosinolates.

For treatment, the first PC contribute more variation in ASCA than in PCA. Also, ASCA gives better clustering.



### Q8: Before doing ASCA with interactions, do you expect a strong interaction effect?

Yes! According to the box plot, there is a strong correlation between time and treatment, especially for PRO, RAPH, GNL, X4OH, GBN, X4MeOH.

### How much of the variation is explained by the interactions?

7.12%

### Answer the question we posed at the start of the exercises: “We would like to know how the different treatments affect the glucosinolate composition and whether the glucosinolate composition varies with time”.

The shoot treatment has a huge positive impact on GBC and NEO, and the root treatment has huge impact on GBN, PRO, GNA and relatively smaller impact on GNL and ALY.

For most of the glucosinolate, the time has a positive impact, exept RAPH and NAS.



*The glucosinolate composition changes due to the different treatments and also due to time. However, the time effect is different for the different treatments.*   听君一席话，胜听一席话



### Q9: This is a point for discussion. What are we looking at when we combine two effects as we did above?

We want to explore the impact one effect on another.

*Here we see the effect of treatment combined with the (time:treatment) interaction. We see that  the JA treatment on ROOT after 14 days returns to Control level (we could interpret this as the plant has recoved, while the SHOOT treatment makes the plant become more different from  control.*

### Make plots of only the treatment effect and only the interaction effect in a similar way as presented above. Discuss the difference with the plot above.

For the plot with only the treatment effect, it only interprets the correlation between treatment groups, without the impact of time.

For the plot with only the interaction, I don't know...



*The combination is obtained by adding the effect of the main treatment factor and the time treatment interaction. The first plot (treatment) shows that the main treatment factor levels. These treatment factor levels are the same for each time point. If in the real data, the treatment effect depends on the time point, then an interaction is able to model this. The interaction shows how the treatment levels change for the different time points. However, the interaction is calculated as the effect on top of the main treatment effect (and therefore the interaction is centered around 0). Combining them gives a better interpretable figure (as seen above). Note that the Y-axis is flipped. If you multiply the PC1 levels with -1, then you can clearly see that the combination is indeed a combination of the main treatment effect and the interaction effect.*
