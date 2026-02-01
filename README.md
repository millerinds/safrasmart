# SafraSmart NDVI (Flask)

Aplicação web para visualizar NDVI por talhão a partir de shapefiles (.zip), com mapa interativo e métricas básicas por polígono.

## Requisitos

- Python 3.9+
- Google Earth Engine (EE) configurado
- Dependências do `requirements.txt`

## Instalação

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Variáveis de ambiente

- `EE_SERVICE_ACCOUNT_JSON` (opcional): JSON completo da service account
- `EE_SERVICE_ACCOUNT_FILE` (opcional): caminho para arquivo JSON da service account
- `PORT` (opcional): porta do Flask (default 8000)
- `CEM_OVERLAY_SIMPLIFY_TOLERANCE` (opcional): tolerância de simplificação do overlay
- `CEM_CLIP_BATCH_THRESHOLD` (opcional): limite de features para recorte em lotes

> Se nenhuma credencial for informada, o app tenta `ee.Initialize()` padrão.

## Como rodar

```bash
python app_flask.py
```

Acesse: `http://localhost:8000`

## Estrutura

- `app_flask.py`: backend Flask + rotas
- `templates/`: HTML do dashboard
- `static/`: CSS
- `smartgreen-dashboard.html`: build estático (referência)

## Observações

- O NDVI é calculado com imagens Sentinel‑2 (L2A) via Earth Engine.
- Valores negativos aparecem em cinza (faixa “sem dado/negativo”).
- O overlay do shapefile é simplificado apenas para performance visual.
