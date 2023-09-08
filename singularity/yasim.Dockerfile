FROM mambaorg/micromamba:1.5.0-alpine
COPY --chown=$MAMBA_USER:$MAMBA_USER yasim-latest.d/env/yasim_minimal.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER yasim-latest.d /yasim-latest.d
COPY --chown=$MAMBA_USER:$MAMBA_USER labw_utils-latest.d /labw_utils-latest.d
COPY --chown=$MAMBA_USER:$MAMBA_USER --chmod=755 pbsim_builder/pbsim /bin/
COPY --chown=$MAMBA_USER:$MAMBA_USER --chmod=755 pbsim_builder/pbsim2 /bin/
COPY --chown=$MAMBA_USER:$MAMBA_USER --chmod=755 pbsim_builder/pbsim3 /bin/

RUN micromamba run -n base python -m pip install --no-deps /labw_utils-latest.d
RUN micromamba run -n base python -m pip install --no-deps /yasim-latest.d

# Check whether installation is successful and versions of installed dependencies.
RUN micromamba run -n base python -m yasim self_check
RUN rm -fr /labw_utils-latest.d /yasim-latest.d
