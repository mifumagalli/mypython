"""
*Use of AI disclaimer*

This code has been developed starting from a template generated with Claude using a prompt asking for a binary classifier of 3D datacubes with masks

"""


import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from astropy.io.ascii import masked
from torch.utils.data import Dataset, DataLoader
import glob
import argparse
from astropy.table import Table
from astropy.io import fits

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from torch.utils.data import Dataset, DataLoader


# Custom 3D Dataset
class ThreeDDataset(Dataset):
    def __init__(self, data, labels):
        """
        Initialize the dataset with 3D data and binary labels

        Args:
        - data: numpy array of shape (num_samples, 1, depth, height, width)
        - labels: numpy array of binary labels (0 or 1)
        """

        self.data = torch.FloatTensor(data)
        self.labels = torch.FloatTensor(labels)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return self.data[idx], self.labels[idx]


# 3D Convolutional Neural Network
class ThreeDClassificationNet(nn.Module):
    def __init__(self, input_shape):
        """
        3D CNN for binary classification

        Args:
        - input_shape: tuple of (depth, height, width)
        """
        super(ThreeDClassificationNet, self).__init__()

        # Convolutional layers
        self.conv_layers = nn.Sequential(
            nn.Conv3d(1, 32, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool3d(2),

            nn.Conv3d(32, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool3d(2),

            nn.Conv3d(64, 128, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool3d(2)
        )

        # Fully connected layers
        self.fc_layers = nn.Sequential(
            nn.Flatten(),
            nn.Linear(128 * (input_shape[0] // 8) * (input_shape[1] // 8) * (input_shape[2] // 8), 256),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(256, 1),
            nn.Sigmoid()  # For binary classification
        )

    def forward(self, x):
        """
        Forward pass through the network

        Args:
        - x: input tensor of shape (batch_size, 1, depth, height, width)
        """
        x = self.conv_layers(x)
        return self.fc_layers(x)


def train_model(model, train_loader, criterion, optimizer, device, epochs=50):
    """
    Train the 3D classification model

    Args:
    - model: PyTorch neural network model
    - train_loader: DataLoader for training data
    - criterion: loss function
    - optimizer: optimization algorithm
    - device: computational device (CPU/GPU)
    - epochs: number of training epochs
    """
    model.train()

    for epoch in range(epochs):
        total_loss = 0

        for batch_data, batch_labels in train_loader:
            # Move data to device
            batch_data = batch_data.to(device)
            batch_labels = batch_labels.unsqueeze(1).to(device)

            # Zero gradients
            optimizer.zero_grad()

            # Forward pass
            outputs = model(batch_data)

            # Compute loss
            loss = criterion(outputs, batch_labels)

            # Backward pass and optimization
            loss.backward()
            optimizer.step()

            total_loss += loss.item()

        # Print average loss per epoch
        print(f'Epoch [{epoch + 1}/{epochs}], Loss: {total_loss / len(train_loader):.4f}')


def evaluate_model(model, test_loader, device):
    """
    Evaluate the model's performance

    Args:
    - model: trained PyTorch neural network
    - test_loader: DataLoader for test data
    - device: computational device

    Returns:
    - Accuracy of the model
    """
    model.eval()
    correct = 0
    total = 0

    with torch.no_grad():
        for batch_data, batch_labels in test_loader:
            batch_data = batch_data.to(device)
            batch_labels = batch_labels.to(device)

            outputs = model(batch_data)
            predicted = (outputs > 0.5).float()

            total += batch_labels.size(0)
            correct += (predicted == batch_labels.unsqueeze(1)).sum().item()

    return correct / total

def cubeloader(pathtodata,batch_size=0):

    """
    Load 3D data and labels from files contained in a directory

    pathtodata: path to directory containing 3D data files

    returns:
        X, y: numpy arrays of 3D data and labels

    """

    # Load data files
    data_files = np.array(glob.glob(pathtodata + '/*.fits'))
    num_samples = len(data_files)

    # Select a subset of files if batch_size > 0
    if(batch_size > 0):
        random_index = np.random.choice(np.arange(num_samples), size=batch_size, replace=False)
        data_files = data_files[random_index]

    X = []
    y = []
    masks=[]

    for file in data_files:
        #read the file
        data = fits.open(file)[0].data

        #need to regularize the data (nan) and make mask
        #thismask = np.isnan(data)
        data[np.isnan(data)] = 0

        #put data between 0 and 1
        data = (data - np.min(data)) / (np.max(data) - np.min(data))
        X.append(data)
        #masks.append(thismask)

        # Extract label from filename
        table = Table.read(file, hdu=1)
        if(('LAE' in table['type']) | ('lae' in table['type'])):
            label=1
        else:
            label=0
        y.append(label)

    #reformat nsample,1,3d_cube_shape
    X = np.array(X)[:,np.newaxis,:,:,:]
    #masks=np.array(masks)[:,np.newaxis,:,:,:]
    #masks = masks.astype(float)
    y = np.array(y)

    return X, y, None



def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Input arguments')
    parser.add_argument('input_path', type=str, help='Path to the data files')
    parser.add_argument('--batch_size', type=int, help='Number of files to select', default=0)
    parser.add_argument('--random_state', type=int, help='Random state for reproducibility', default=42)

    args = parser.parse_args()

    # Set random seed for reproducibility
    torch.manual_seed(args.random_state)
    np.random.seed(args.random_state)
    
    # Determine device
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    print(f'Using device: {device}')

    #load data
    X, y, masks= cubeloader(args.input_path, batch_size=args.batch_size)
    num_samples = len(y)

    # Split data into training and testing
    train_ratio = 0.8
    split_idx = int(num_samples * train_ratio)
    
    X_train, X_test = X[:split_idx], X[split_idx:]
    y_train, y_test = y[:split_idx], y[split_idx:]

    # Create datasets and dataloaders
    train_dataset = ThreeDDataset(X_train, y_train)
    test_dataset = ThreeDDataset(X_test, y_test)

    train_loader = DataLoader(train_dataset, batch_size=16, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=16, shuffle=False)


    # Initialize model
    input_shape=X_train.shape[2:]
    model = ThreeDClassificationNet(input_shape)
    if torch.cuda.device_count() > 1:
        print("Let's use", torch.cuda.device_count(), "GPUs!")
        # dim = 0 [30, xxx] -> [10, ...], [10, ...], [10, ...] on 3 GPUs
        model = nn.DataParallel(model)
    model.to(device)

    # Define loss and optimizer
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    # Train the model
    train_model(model, train_loader, criterion, optimizer, device)
    
    # Evaluate the model
    accuracy = evaluate_model(model, test_loader, device)
    print(f'Test Accuracy: {accuracy * 100:.2f}%')
    
    # Example prediction
    #sample_data = torch.FloatTensor(np.random.randn(1, 1, *input_shape)).to(device)
    #with torch.no_grad():
    #    prediction = model(sample_data)
    #    result = "Yes" if prediction.item() > 0.5 else "No"
    #    print(f'Sample Prediction: {result} (Confidence: {prediction.item():.2f})')

if __name__ == "__main__":
    main()
