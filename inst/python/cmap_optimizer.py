import os,sys
import numpy as np
import pandas as pd
import cv2
import logging
import torch
from torch.nn.functional import softmax, cosine_similarity
sys.path.insert(1,os.path.join(os.path.abspath('.'),"pytorch-msssim"))
import pytorch_msssim_me as pytorch_msssim

def infor(data):
    a = pd.value_counts(data) / len(data)
    return sum(np.log2(a) * a * (-1))

class Mapping:

    def __init__(
        self,
        C,
        S,
        device="cpu",
        random_seed_set=None
    ):
        """
        
        """
        self.C = torch.tensor(C, device=device, dtype=torch.float32)
        self.S = torch.tensor(S, device=device, dtype=torch.float32)

        self.random_seed_set = random_seed_set
        if self.random_seed_set:
            np.random.seed(seed=self.random_seed_set)
        self.M = np.random.normal(0, 1, (C.shape[0], S.shape[0]))
        self.M = torch.tensor(
            self.M, device=device, requires_grad=True, dtype=torch.float32
        )

        self.S_n = (255*(self.S - self.S.min())/(self.S.max() - self.S.min()))
        self.S_n = self.S_n.repeat(1,1,1,1)


    def _loss_fn(self, para_distance, para_density, verbose=True):
        """

        """
        M_probs = softmax(self.M, dim=1)

        S_pred_t = torch.matmul(M_probs.t(), self.C)

        S_pred_n = (255*(S_pred_t - S_pred_t.min())/(S_pred_t.max() - S_pred_t.min()))
        S_pred_n = S_pred_n.repeat(1,1,1,1)

        map_spot_index = torch.max(M_probs, 1)[1].cpu().numpy()
        H_info = infor(map_spot_index)
        total_loss = para_distance*(1 - pytorch_msssim.ssim(S_pred_n, self.S_n, window_size=2)) - para_density*H_info

        return (total_loss)
    


    def train(self, num_epochs, para_distance, para_density, learning_rate=0.1, print_each=100):
        """
        Function: Optimize the mapping matrix to find the global optimization.

        Parameters:
            para_distance: The weight for the SSIM term.
            para_density: The weight for the entropy term.

        Return:
            a cell-by-spot corresponding matrix (AnnData type), containing the probability of mapping cell i to spot j.

        """
        if self.random_seed_set:
            torch.manual_seed(seed=self.random_seed_set)
        optimizer = torch.optim.NAdam([self.M], lr=learning_rate)


        if print_each:
            logging.info(f"Printing scores every {print_each} epochs.")
            

        for t in range(num_epochs):
            if print_each is None or t % print_each != 0:
                run_loss = self._loss_fn(para_distance,para_density,verbose=False)
            else:
                run_loss = self._loss_fn(para_distance,para_density,verbose=True)

            print(run_loss)

            loss = run_loss

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        with torch.no_grad():
            output = softmax(self.M, dim=1).cpu().numpy()
            return output


