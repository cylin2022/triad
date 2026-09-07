# TRIAD Platform — port 5001, v1.1
#
# This image ships WITHOUT netMHCpan-4.2 itself: its academic license from
# DTU Health Tech forbids redistribution to third parties. Each user must
# obtain their own copy from https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/
# and install it after starting the container — either by mounting a volume
# at NETMHC_HOME, or through the built-in /setup upload wizard.
#
# The netMHCpan-4.2 binary is a statically-linked ELF executable, so no
# tcsh/awk runtime is required here (the platform calls it directly).

FROM python:3.10-slim

LABEL org.opencontainers.image.title="triad-platform" \
      org.opencontainers.image.version="v1.1"

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .
RUN chmod +x entrypoint.sh

# Mount point for the user-supplied netMHCpan-4.2 package.
ENV NETMHC_HOME=/opt/netMHCpan-4.2
RUN mkdir -p "$NETMHC_HOME" uploads results logs

ENV FLASK_PORT=5001
EXPOSE 5001

ENTRYPOINT ["./entrypoint.sh"]
CMD ["python3", "app.py"]
