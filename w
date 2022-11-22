{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "c2b8caf5",
   "metadata": {},
   "outputs": [],
   "source": [
    "import wandb\n",
    "api = wandb.Api()\n",
    "\n",
    "# run is specified by <entity>/<project>/<run id>\n",
    "run = api.run(\"birds/yeast_pvalue_otherstats/8blfajvj\")\n",
    "run1= api.run(\"birds/yeast_pvalue_otherstats/1hrpwn4b\")\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "67944278",
   "metadata": {},
   "outputs": [],
   "source": [
    "# save the metrics for the run to a csv file\n",
    "metrics_dataframe = run.history()\n",
    "metrics_dataframe.to_csv(\"groundtruth0.csv\")\n",
    "metrics_dataframe1 = run1.history()\n",
    "metrics_dataframe1.to_csv(\"groundtruth1.csv\")"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.7"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
