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
        - data: numpy array of shape (num_samples, depth, height, width)
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
            nn.Linear(128 * (input_shape[0]//8) * (input_shape[1]//8) * (input_shape[2]//8), 256),
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
        print(f'Epoch [{epoch+1}/{epochs}], Loss: {total_loss/len(train_loader):.4f}')

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

def main():
    # Set random seed for reproducibility
    torch.manual_seed(42)
    np.random.seed(42)
    
    # Determine device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Using device: {device}')
    
    # Generate sample 3D data (replace with your actual data)
    num_samples = 1000
    input_shape = (16, 32, 32)  # depth, height, width
    
    # Simulate 3D data with random noise and labels
    X = np.random.randn(num_samples, 1, *input_shape)
    y = np.random.randint(2, size=num_samples).astype(float)
    
    # Split data into training and testing
    train_ratio = 0.8
    split_idx = int(num_samples * train_ratio)
    
    X_train, X_test = X[:split_idx], X[split_idx:]
    y_train, y_test = y[:split_idx], y[split_idx:]
    
    # Create datasets and dataloaders
    train_dataset = ThreeDDataset(X_train, y_train)
    test_dataset = ThreeDDataset(X_test, y_test)
    
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
    
    # Initialize model
    model = ThreeDClassificationNet(input_shape).to(device)
    
    # Define loss and optimizer
    criterion = nn.BCELoss()  # Binary Cross Entropy for binary classification
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    # Train the model
    train_model(model, train_loader, criterion, optimizer, device)
    
    # Evaluate the model
    accuracy = evaluate_model(model, test_loader, device)
    print(f'Test Accuracy: {accuracy * 100:.2f}%')
    
    # Example prediction
    sample_data = torch.FloatTensor(np.random.randn(1, 1, *input_shape)).to(device)
    with torch.no_grad():
        prediction = model(sample_data)
        result = "Yes" if prediction.item() > 0.5 else "No"
        print(f'Sample Prediction: {result} (Confidence: {prediction.item():.2f})')

if __name__ == "__main__":
    main()
