# Tiny Bit of Background

Here we calculate two-sample bounds based on simulated data. This allows us to explore the behavior of two-sample bounds when multiple instruments influence the exposure. The model used is relatively flexible and allows us to include both valid and invalid instruments, and a pair of dependent instruments. 


# Simulation Scheme

The data is simulated in the following way:

* A value of the confounder $U$ is drawn from a standard Normal distribution

* Values of independent IVs ($Z_1, ..., Z_p$) are drawn based on user specified probabilities.

* Optionally, a pair of dependent IVs, $Z_{D1}, Z_{D2}$, are simulated as follows:
    * User specifies $\rho$
    * $(N_{D1}, N_{D2}) \sim N\left(\begin{pmatrix} 0 \\ 0 \end{pmatrix},\begin{pmatrix} 1 & \rho \\ \rho & 1 \end{pmatrix} \right)$, then $$Z_{Di} = \left\{\begin{array}{cl} 0 & \text{ if } N_{Di} \le z_{p_1} \\ 1 & \text{ if } z_{p_0} < N_{Di} \le z_{p_0 + p_1} \\ 2 & \text{ otherwise} \end{array}  \right. ,$$
      where $z_k$ is the $k'th$ percentile of the standard normal ($P(Z \le z_k) = k)$ when $Z \sim N(0,1)$). 

* A binary treatment is simulated as $X | Z_1,...,Z_p, Z_{D1}, Z_{D2}, U \sim \text{Bernoulli}(p_X)$, where $$\text{logit}(p_X) = \gamma_0 + \gamma_1 \cdot Z_1 + ... + \gamma_p \cdot Z_p + \gamma_{D1} \cdot Z_{D1} + \gamma_{D2} \cdot Z_{D2} + \gamma_U \cdot U.$$ The parameters $\gamma_1, ..., \gamma_p, \gamma_U, \gamma_{D1}, \gamma_{D2}$ are all user specified, while $\gamma$ is set to $-\sum_{j = 1}^p \gamma_j$.

* A binary outcome is simulated as $Y | X, U, Z_1,...,Z_p, Z_{D1}, Z_{D2}, U \sim \text{Bernoulli}(p_Y)$, where $$\text{logit}(p_Y) = \beta_0 + \beta_X \cdot X + \beta_U \cdot U + \beta_1 \cdot Z_1 + ... + \beta_p \cdot Z_p + \beta_{D1} \cdot Z_{D1} + \beta_{D2} \cdot Z_{D2}.$$ The parameters $\beta_X, \beta_U, \beta_1, ..., \beta_p, \beta_{D1}, \beta_{D2}$ are all user specified, while $\beta_0$ is set to $-\beta_X/2$. When no invalid IVs are included, $\beta_1, ..., \beta_p, \beta_{D1}, \beta_{D2}$ are all $0$.

* The above steps are repeated to create a data set of the sample size specified. 

# Results

Based on the simulated data, values of $P(X = 1 | Z_j = z)$ and $P(Y = 1 | Z_j = z)$ are calculated and used to find two-sample bounds. These bounds are found using the function [`ATEBounds::get_bounds()`](https://rmtrane.github.io/ATEBounds/reference/get_bounds.html). You will find the bounds on a figure under the tab "Figure", in a table under the tab "Table", and the simulated data under the tab "Simulated Data". 
