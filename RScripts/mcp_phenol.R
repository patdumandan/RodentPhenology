# data
dat_list=list(N=length(PB_gam_con$X),
              proportion=PB_gam_con$proportion,
              month=PB_gam_con$month,
              year=PB_gam_con$years)

#mcp

model = list(proportion
                        ~1,
                        ~0+month,
                        ~1+month
            )  # three intercept-only segments
fit_mcp = mcp(model, data = PB_gam_con, par_x = "month")

summary(fit_mcp)
plot(fit_mcp) + plot_pars(fit_mcp, pars = c("cp_1","cp_2"), type = "dens_overlay")
