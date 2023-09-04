FROM mambaorg/micromamba:alpine-1.5.0
COPY --chown=$MAMBA_USER:$MAMBA_USER yasim-latest.d/env/yasim_minimal.yml /tmp/env.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER --chmod=755 pbsim_builder/pbsim /bin/
COPY --chown=$MAMBA_USER:$MAMBA_USER --chmod=755 pbsim_builder/pbsim2 /bin/
COPY --chown=$MAMBA_USER:$MAMBA_USER --chmod=755 pbsim_builder/pbsim3 /bin/
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes \
RUN python yasim-latest.d/setup.py install
RUN python -m yasim self_check
