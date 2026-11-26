"""Catalogue of ten deep-learning architectures."""
DEEP_MODELS=("mlp","residual_mlp","cnn1d","lstm","gru","transformer","autoencoder","variational_autoencoder","denoising_autoencoder","multi_task_mlp")
from .neural_network import build_model,train_torch_model
