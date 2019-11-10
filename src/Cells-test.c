/*************************************************************************
  Unit tests for routines in Cells.h 
*************************************************************************/

#include "Cells.h"

#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define RESET   "\x1b[0m"

int checkn_Nbrs();

int CheckNodeType();

unsigned long CompareLinks(enum nodeType_t type);

int main() {
    unsigned long int nodeIt;

    char *data = "Networks/10x10.dat";
    char *network = "Networks/10x10.net";
    char *coords = "Networks/10x10.cor";
    int return_Val = 0;
    unsigned long int *sizes = NULL;
    printf("Read in & Verify 10^2 grid     : ");
    return_Val = ReadLattice(data, network, coords);
    if (return_Val == 0)
        return_Val = VerifyState();
    if (return_Val == 0) {
        printf(GREEN "PASSED\n" RESET);
    } else {
        printf(RED "FAILED\n" RESET);
    }

    if (return_Val == 0) {
        printf("Check n_Nbrs : ");
        if (checkn_Nbrs() == 0)
            printf(GREEN "PASSED\n" RESET);
        else
            printf(RED "FAILED\n" RESET);
        fflush(stdout);
        struct node_t *center = lttc->nodes[lttc->center_Idx];
        printf("Verifying 0 Tumors with GetNetworkDetails() :");
        sizes = GetNetworkDetails(lttc->cells[CANCER], lttc->n_Cells[CANCER]);
        if (sizes[0] == 0)
            printf(GREEN "PASSED\n" RESET);
        else
            printf(RED "FAILED\n" RESET);
        free(sizes);


        fflush(stdout);

        printf("Seeding 2 Tumors & Verify Count: ");
        SeedCancer(lttc->nodes[0], 2, 2);
        SeedCancer(center, 3, 2);
        if (lttc->n_Cells[CANCER] == 26 &&
            lttc->n_Cells[NORMAL] == 74 &&
            lttc->n_Cells[INFECT_INT] == 0 &&
            lttc->n_Cells[INFECT_EXT] == 0)
            printf(GREEN "PASSED\n" RESET);
        else {
            printf(RED "FAILED\n" RESET);
            printf("------Cancer (26): %ld\n", lttc->n_Cells[CANCER]);
            printf("------Normal (74): %ld\n", lttc->n_Cells[NORMAL]);
            printf("------INFECT_INT (0): %ld\n", lttc->n_Cells[INFECT_INT]);
            printf("------INFECT_EXT (0): %ld\n", lttc->n_Cells[INFECT_EXT]);
        }


        fflush(stdout);

        printf("Verifying Tumors with GetNetworkDetails() :");
        sizes = GetNetworkDetails(lttc->cells[CANCER], lttc->n_Cells[CANCER]);
        if ((sizes[0] == 2) && (sizes[1] == 13) && (sizes[2] == 13))
            printf(GREEN "PASSED\n" RESET);
        else
            printf(RED "FAILED\n" RESET);
        free(sizes);
        fflush(stdout);

        printf("Verifying n_Targets[INFECT] :");
        if (lttc->n_Target_Nodes[INFECT_INT][0] == 0 &&
            lttc->n_Target_Nodes[INFECT_INT][1] == 0 &&
            lttc->n_Target_Nodes[INFECT_INT][2] == 0 &&
            lttc->n_Target_Nodes[INFECT_INT][3] == 0 &&
            lttc->n_Target_Nodes[INFECT_INT][4] == 0 &&
            lttc->n_Target_Nodes[INFECT_EXT][0] == 0 &&
            lttc->n_Target_Nodes[INFECT_EXT][1] == 0 &&
            lttc->n_Target_Nodes[INFECT_EXT][2] == 0 &&
            lttc->n_Target_Nodes[INFECT_EXT][3] == 0 &&
            lttc->n_Target_Nodes[INFECT_EXT][4] == 0)
            printf(GREEN "PASSED\n" RESET);
        else
            printf(RED "FAILED\n" RESET);
        fflush(stdout);

        printf("Infect 2 Tumors & Verify Count: ");
        VirusBallAroundNode(center, 13);
        VirusBallAroundNode(lttc->nodes[0], 13);
        if (lttc->n_Cells[INFECT_INT] + lttc->n_Cells[INFECT_EXT] == 26)
            printf(GREEN "PASSED\n" RESET);
        else
            printf(RED "FAILED\n" RESET);
        fflush(stdout);

        printf("Verifying with GetNetworkDetails() :");
        sizes = GetNetworkDetails(lttc->cells[CANCER], lttc->n_Cells[CANCER]);
        if (sizes[0] == 0)
            printf(GREEN "PASSED\n" RESET);
        else
            printf(RED "FAILED\n" RESET);
        free(sizes);
        fflush(stdout);

        printf("Verifying with n_Targets[INFECT] :");
        if (lttc->n_Target_Nodes[INFECT_INT][0] == 0 &&
            lttc->n_Target_Nodes[INFECT_INT][1] == 0 &&
            lttc->n_Target_Nodes[INFECT_INT][2] == 0 &&
            lttc->n_Target_Nodes[INFECT_INT][3] == 0 &&
            lttc->n_Target_Nodes[INFECT_INT][4] == 0 &&
            lttc->n_Target_Nodes[INFECT_EXT][0] == 0 &&
            lttc->n_Target_Nodes[INFECT_EXT][1] == 0 &&
            lttc->n_Target_Nodes[INFECT_EXT][2] == 0 &&
            lttc->n_Target_Nodes[INFECT_EXT][3] == 0 &&
            lttc->n_Target_Nodes[INFECT_EXT][4] == 0)
            printf(GREEN "PASSED\n" RESET);
        else
            printf(RED "FAILED\n" RESET);
        fflush(stdout);

        printf("Seed 2x13node tumors, then infect 14 nodes.  Verify Count & n_Nbrs:");
        SeedCancer(lttc->nodes[0], 2, 2);
        SeedCancer(center, 3, 2);
        VirusBallAroundNode(center, 14);
        if (lttc->n_Cells[INFECT_INT] + lttc->n_Cells[INFECT_EXT] == 14 &&
            checkn_Nbrs() == 0)
            printf(GREEN "PASSED\n" RESET);
        else
            printf(RED "FAILED\n" RESET);
        fflush(stdout);


        printf("Testing Viral Growth from 1 node to full & verify n_Links & target types: ");
        SeedCancer(center, 3, 10 * 10);
        VirusBallAroundNode(center, 1);
        rates[NORMAL][GROW] = 0;
        rates[NORMAL][KILL] = 0;
        rates[CANCER][GROW] = 0;
        rates[CANCER][KILL] = 0;
        rates[INFECT_INT][GROW] = 10;
        rates[INFECT_INT][KILL] = 0;
        rates[INFECT_EXT][GROW] = 10;
        rates[INFECT_EXT][KILL] = 0;
        unsigned long It;
        for (It = 1; It < 10 * 10; It++) {
            Simulate();

            if (CheckNodeType() != 0)
                break;
            if (CompareLinks(INFECT_EXT) + CompareLinks(INFECT_INT) != 0)
                break;
            if (checkn_Nbrs() != 0)
                break;
        }
        if (lttc->n_Cells[INFECT_INT] + lttc->n_Cells[INFECT_EXT] == 10 * 10)
            printf(GREEN "PASSED\n" RESET);
        else
            printf(RED "FAILED\n" RESET);
        fflush(stdout);

        rates[NORMAL][GROW] = 0;
        rates[NORMAL][KILL] = 0;
        rates[CANCER][GROW] = 10;
        rates[CANCER][KILL] = 0;
        rates[INFECT_INT][GROW] = 0;
        rates[INFECT_INT][KILL] = 0;
        rates[INFECT_EXT][GROW] = 0;
        rates[INFECT_EXT][KILL] = 0;
        printf("Testing Cancer Growth Next to Virus: ");
        SeedCancer(center, 10 * 10, 1);
        VirusBallAroundNode(center, 13);
        ChangeNodeType(lttc->nodes[41], VACANT);
        ChangeNodeType(lttc->nodes[32], VACANT);
        ChangeNodeType(lttc->nodes[23], VACANT);
        ChangeNodeType(lttc->nodes[14], VACANT);
        ChangeNodeType(lttc->nodes[25], VACANT);
        ChangeNodeType(lttc->nodes[36], VACANT);
        ChangeNodeType(lttc->nodes[47], VACANT);
        ChangeNodeType(lttc->nodes[56], VACANT);
        ChangeNodeType(lttc->nodes[65], VACANT);
        ChangeNodeType(lttc->nodes[74], VACANT);
        ChangeNodeType(lttc->nodes[63], VACANT);
        ChangeNodeType(lttc->nodes[52], VACANT);
        {
            int error = 0;
            while (lttc->n_Cells[CANCER] < 87) {
                Simulate();
                if (VerifyState() == 0 &&
                    checkn_Nbrs() == 0 &&
                    CheckNodeType() == 0 &&
                    CompareLinks(INFECT_INT) + CompareLinks(INFECT_EXT) == 0)
                    continue;
                else {
                    printf(RED "FAILED\n" RESET);
                    error = 1;
                    break;
                }

            }
            if (error == 0)
                printf(GREEN "PASSED\n" RESET);
            fflush(stdout);
        }

        {
            printf("Changing nodes to Vacant then to Cancer: ");
            int error = 0;
            for (nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++) {
                struct node_t *node = lttc->nodes[nodeIt];
                ChangeNodeType(node, VACANT);
                ChangeNodeType(node, CANCER);
                if (VerifyState() == 0)
                    continue;
                else {
                    printf(RED "FAILED\n" RESET);
                    error = 1;
                    break;
                }
            }
            if (error == 0)
                printf(GREEN "PASSED\n" RESET);
            fflush(stdout);
        }
        {
            printf("Testing Normal Growth next to target: ");
            for (nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++) {
                struct node_t *node = lttc->nodes[nodeIt];
                ChangeNodeType(node, NORMAL);
            }
            ChangeNodeType(lttc->nodes[0], VACANT);
            ChangeNodeType(lttc->nodes[1], VACANT);
            ChangeNodeType(lttc->nodes[1], NORMAL);
            if (VerifyState() == 0)
                printf(GREEN "PASSED\n" RESET);
            else
                printf(RED "FAILED\n" RESET);
            fflush(stdout);
        }
        {

            printf("Random Infection Count:");
            percent_Infected = 20;
            infect_Type = RANDOM;
            for (nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++)
                ChangeNodeType(lttc->nodes[nodeIt], CANCER);
            InjectVirus();
            sizes = GetNetworkDetails(lttc->cells[INFECT_INT], lttc->n_Cells[INFECT_INT]);
            if (sizes[0] > 10 && lttc->n_Cells[INFECT_INT] == (lttc->n_Nodes * (int) percent_Infected) / 100)
                printf(GREEN "PASSED\n" RESET);
            else
                printf(RED "FAILED\n" RESET);
            fflush(stdout);
            free(sizes);

            printf("Center Infection Count:");
            for (nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++)
                ChangeNodeType(lttc->nodes[nodeIt], CANCER);
            infect_Type = CENTER;
            InjectVirus();
            sizes = GetNetworkDetails(lttc->cells[INFECT_INT], lttc->n_Cells[INFECT_INT]);
            if (sizes[0] == 1 && sizes[1] == lttc->n_Nodes * (int) percent_Infected / 100)
                printf(GREEN "PASSED\n" RESET);
            else
                printf(RED "FAILED\n" RESET);
            fflush(stdout);
            free(sizes);

            printf("Multi Infection Count:");
            percent_Infected = (float) 18.75;
            for (nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++) {
                enum nodeType_t type = CANCER;
                if (lttc->nodes[nodeIt]->pos->data[0] == 1 ||
                    lttc->nodes[nodeIt]->pos->data[0] == 10 ||
                    lttc->nodes[nodeIt]->pos->data[1] == 1 ||
                    lttc->nodes[nodeIt]->pos->data[1] == 10)
                    type = NORMAL;

                ChangeNodeType(lttc->nodes[nodeIt], type);
            }
            infect_Type = MULTINODE;
            InjectVirus();
            sizes = GetNetworkDetails(lttc->cells[INFECT_INT], lttc->n_Cells[INFECT_INT]);
            int passed = 0;
            for (It = 1; It <= sizes[0]; It++) {
                if (sizes[It] % (int) percent_Infected / 3 / 100 * 64 != 0) {
                    passed = 0;
                    break;
                } else
                    passed = 1;
            }
            free(sizes);

            unsigned long n_Infected = lttc->n_Cells[INFECT_EXT] +
                                       lttc->n_Cells[INFECT_INT];
            if (passed == 1 && n_Infected == (int) (percent_Infected / 100 * 64))
                printf(GREEN "PASSED\n" RESET);
            else
                printf(RED "FAILED\n" RESET);
            fflush(stdout);
            printf("Perimeter Infection Count:");
            percent_Infected = (float) 18.75;
            for (nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++) {
                enum nodeType_t type = CANCER;
                if (lttc->nodes[nodeIt]->pos->data[0] == 1 ||
                    lttc->nodes[nodeIt]->pos->data[0] == 10 ||
                    lttc->nodes[nodeIt]->pos->data[1] == 1 ||
                    lttc->nodes[nodeIt]->pos->data[1] == 10)
                    type = NORMAL;

                ChangeNodeType(lttc->nodes[nodeIt], type);
            }
            infect_Type = PERIMETER;
            InjectVirus();
            n_Infected = lttc->n_Cells[INFECT_EXT] +
                         lttc->n_Cells[INFECT_INT];
            if (n_Infected == (int) (percent_Infected / 100 * 64))
                printf(GREEN "PASSED\n" RESET);
            else
                printf(RED "FAILED\n" RESET);
            fflush(stdout);

        }
        DeleteLattice(lttc);
    }

    /* Read in 10^3 Grid */
    data = "Networks/10x10x10.dat";
    network = "Networks/10x10x10.net";
    coords = "Networks/10x10x10.cor";

    printf("Read in & Verify 10^3 grid     : ");
    return_Val = ReadLattice(data, network, coords);
    if (return_Val == 0)
        return_Val = VerifyState();
    if (return_Val == 0) {
        printf(GREEN "PASSED\n" RESET);
    } else {
        printf(RED "FAILED\n" RESET);
        return 0;
    }
    printf("Perimeter Infection Count:");
    percent_Infected = 104.0 / 512.0 * 100.0;
    for (nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++) {
        enum nodeType_t type = CANCER;
        if (lttc->nodes[nodeIt]->pos->data[0] == 1 ||
            lttc->nodes[nodeIt]->pos->data[0] == 10 ||
            lttc->nodes[nodeIt]->pos->data[1] == 1 ||
            lttc->nodes[nodeIt]->pos->data[1] == 10 ||
            lttc->nodes[nodeIt]->pos->data[2] == 1 ||
            lttc->nodes[nodeIt]->pos->data[2] == 10)
            type = NORMAL;

        ChangeNodeType(lttc->nodes[nodeIt], type);
    }
    infect_Type = PERIMETER;
    InjectVirus();

    unsigned long n_Infected = lttc->n_Cells[INFECT_EXT] +
                               lttc->n_Cells[INFECT_INT];
    if (n_Infected == (int) (percent_Infected / 100 * 512))
        printf(GREEN "PASSED\n" RESET);
    else
        printf(RED "FAILED\n" RESET);
    fflush(stdout);
    OutputLatticeState("output.csv");
    printf("Multi Infection Count:");
    /* This also tests correct behaviour of VirusBallAroundNode( node, N )
       when node is already infected. */
    percent_Infected = (float) 18.75;
    for (nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++) {
        enum nodeType_t type = CANCER;
        if (lttc->nodes[nodeIt]->pos->data[0] == 1 ||
            lttc->nodes[nodeIt]->pos->data[0] == 10 ||
            lttc->nodes[nodeIt]->pos->data[1] == 1 ||
            lttc->nodes[nodeIt]->pos->data[1] == 10 ||
            lttc->nodes[nodeIt]->pos->data[2] == 1 ||
            lttc->nodes[nodeIt]->pos->data[2] == 10)
            type = NORMAL;

        ChangeNodeType(lttc->nodes[nodeIt], type);
    }
    infect_Type = MULTINODE;
    InjectVirus();
    sizes = GetNetworkDetails(lttc->cells[INFECT_EXT], lttc->n_Cells[INFECT_EXT]);

    int passed = 0;
    unsigned long It;
    for (It = 1; It <= sizes[0]; It++) {
        printf("size of %ld on %ld\n", sizes[It], It);
        if (sizes[It] % (int) percent_Infected / 3 / 100 * 512 != 0) {
            passed = 0;
            break;
        } else
            passed = 1;
    }
    free(sizes);

    n_Infected = lttc->n_Cells[INFECT_EXT] +
                 lttc->n_Cells[INFECT_INT];
    if (passed == 1 && n_Infected + 1 == (int) (percent_Infected / 100 * 512))
        /*  100+(10+10+16)*8+100 = 488 NORMAL Cells
            1000-488     = 512 CANCER Cells */
        printf(GREEN "PASSED\n" RESET);
    else {
        printf(RED "FAILED\n" RESET);
        printf("------Passed: %d\n", passed);
        printf("------N_Infect: %lu\n", n_Infected);
        printf("------percent_Infect: %f\n", percent_Infected);
        printf("------512: %d\n", (int) (percent_Infected / 100 * 512));
    }
    fflush(stdout);

    printf("Multi Infection Check Center:");
    Point *point = CenterOfMass(lttc->cells[INFECT_EXT], lttc->n_Cells[INFECT_EXT]);
    passed = 0;
    for (It = 0; It < point->N; It++) {
        if (fabs(5.5 - point->data[It]) > 1.5) {
            passed = 0;
            break;
        } else passed = 1;
    }
    if (passed == 1)
        printf(GREEN "PASSED\n" RESET);
    else
        printf(RED "FAILED\n" RESET);
    fflush(stdout);

    DeletePoint(point);
    DeleteLattice(lttc);

    data = "Networks/3d_grid_10.dat";
    network = "Networks/3d_grid_10.net";
    coords = "Networks/3d_grid_10.state";

    printf("Read in & Verify 10^3 grid     : ");
    return_Val = ReadLattice(data, network, coords);
    if (return_Val == 0)
        return_Val = VerifyState();
    if (return_Val == 0) {

        printf(GREEN "PASSED\n" RESET);
    } else {
        printf(RED "FAILED\n" RESET);
    }
    DeleteLattice(lttc);

    return 0;
}

int checkn_Nbrs() {
    int typeIt, nbrIt;
    unsigned long int nodeIt;
    struct node_t *node;
    for (nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++) {
        int n_Nbrs[N_CELL_TYPES + 1] = {0, 0, 0, 0};
        node = lttc->nodes[nodeIt];
        for (nbrIt = 0; nbrIt < node->n_Neighbors; nbrIt++)
            n_Nbrs[node->neighbors[nbrIt]->type]++;

        for (typeIt = 0; typeIt < N_CELL_TYPES; typeIt++) {
            if (n_Nbrs[typeIt] != node->n_Nbrs[typeIt]) {
                fprintf(stderr, "n_Nbrs count is off\n");
                return (-1);
            }
        }
    }

    return 0;
}

int CheckNodeType() {
    int typeIt;
    for (typeIt = NORMAL; typeIt < N_CELL_TYPES; typeIt++) {
        unsigned long n_NbrsIt;
        for (n_NbrsIt = 0; n_NbrsIt <= lttc->n_Nbrs_Max; n_NbrsIt++) {
            unsigned long int nodeIt;
            for (nodeIt = 0; nodeIt < lttc->n_Target_Nodes[typeIt][n_NbrsIt]; nodeIt++) {
                if (lttc->target_Nodes[typeIt][n_NbrsIt][nodeIt]->type != TargetOf(typeIt)) {
                    printf("\nWrong node type in target_Nodes; count: %lu\n", lttc->count);
                    return (-1);
                }
            }
        }
    }
    return 0;
}

unsigned long CompareLinks(enum nodeType_t type) {
    unsigned long int IC_Links_a = 0;
    unsigned long int nodeIt;
    /* Count n_links from lttc->cells[] */
    for (nodeIt = 0; nodeIt < lttc->n_Cells[TargetOf(type)]; nodeIt++) {
        unsigned long int nbrIt;
        for (nbrIt = 0; nbrIt < lttc->cells[TargetOf(type)][nodeIt]->n_Neighbors; nbrIt++) {
            if (lttc->cells[TargetOf(type)][nodeIt]->neighbors[nbrIt]->type == type)
                IC_Links_a++;
        }
    }
    unsigned long int IC_Links_b = 0;
    unsigned long int n_nbrIt;
    /* Count n_links from n_Target_Nodes[] */
    for (n_nbrIt = 0; n_nbrIt <= lttc->n_Nbrs_Max; n_nbrIt++) {
        IC_Links_b += lttc->n_Target_Nodes[type][n_nbrIt] * n_nbrIt;
    }
    /* Return Difference */
    return IC_Links_a - IC_Links_b;
}


/** VerifyState()******************************************
* Check for some inconsistancies in the lattice state.
* Useful to check changes in code haven't cause problems
**************************************************************/
int VerifyState() {
    struct node_t *node;
    int n_Nbrs[N_CELL_TYPES + 1] = {0, 0, 0, 0};
    unsigned long int n_Cells[N_CELL_TYPES + 1] = {0, 0, 0, 0};
    unsigned long int *n_Target_Nodes[N_CELL_TYPES + 1];
    unsigned long int nbrIt;
    int typeIt;
    unsigned long int nodeIt;
    unsigned long int count = 0;

    for (typeIt = 0; typeIt < N_CELL_TYPES; typeIt++) {
        n_Target_Nodes[typeIt] = calloc((lttc->n_Nbrs_Max + 1), sizeof(**n_Target_Nodes));
    }

    for (nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++) {
        node = lttc->nodes[nodeIt];
        count++;
        /* Assume node->type,
  n_Neighbors,
  neighbors */

        /* Verify node->n_Nbrs[] is correct */
        for (nbrIt = 0; nbrIt < node->n_Neighbors; nbrIt++)
            n_Nbrs[node->neighbors[nbrIt]->type]++;
        for (typeIt = 0; typeIt < N_CELL_TYPES; typeIt++) {
            if (n_Nbrs[typeIt] != node->n_Nbrs[typeIt]) {
                fprintf(stderr, "n_Nbrs count is off\n");
                return (-1);
            }
        }

        /* Verify list_Idx is correct */
        if (node->list_Idx != NO_IDX &&
            lttc->cells[node->type][node->list_Idx] != node) {
            fprintf(stderr, "Index to cells[] is incorrect\n");
            return (-1);
        }

        /* Verify targets_Idx[] are correct */
        for (typeIt = 0; typeIt < N_CELL_TYPES; typeIt++) {
            if (node->targets_Idx[typeIt] == NO_IDX) {
            } else {
                /* Verify correct address at Idx */
                if (lttc->target_Nodes[typeIt]
                    [n_Nbrs[typeIt]]
                    [node->targets_Idx[typeIt]] != node) {
                    fprintf(stderr, "Index to target_Nodes[] is incorrect\n");
                    return (-1);
                } else {
                    if (n_Nbrs[typeIt] == 0) {
                        fprintf(stderr, "Node w/ nonzero target_Idx has no neighbors that \
can replicate to it.\n");
                        return (-1);
                    }

                    n_Target_Nodes[typeIt][n_Nbrs[typeIt]]++;
                }
            }
        }

        n_Cells[node->type]++;
        n_Nbrs[0] = 0;
        n_Nbrs[1] = 0;
        n_Nbrs[2] = 0;
    }

    /* Verify lttc->n_Cells[] is correct */
    for (typeIt = 0; typeIt < N_CELL_TYPES; typeIt++) {
        if (n_Cells[typeIt] != lttc->n_Cells[typeIt]) {
            fprintf(stderr, "n_Cells[] is incorrect\n");
            return (-1);
        }
    }

    /* Verify cells[] is correct */
    for (typeIt = 0; typeIt < N_CELL_TYPES; typeIt++) {
        for (nodeIt = 0; nodeIt < n_Cells[typeIt]; nodeIt++) {
            if (lttc->cells[typeIt][nodeIt]->type != typeIt) {
                fprintf(stderr, "Cells of wrong type in cells[].\n");
                return (-1);
            }
        }
        /* Make sure all elements after n_Cells[typeIt] are NULL */
        for (nodeIt = n_Cells[typeIt]; nodeIt < lttc->n_Nodes; nodeIt++) {
            if (lttc->cells[typeIt][nodeIt] != NULL) {
                fprintf(stderr, "Not NULL cells exist in cells[%d] after n_Cells.\n",
                        typeIt);
                return (-1);
            }
        }
    }

    /* Verify lttc->n_Target_Nodes is correct */
    for (typeIt = 0; typeIt < N_CELL_TYPES; typeIt++) {
        for (nbrIt = 0; nbrIt < lttc->n_Nbrs_Max; nbrIt++) {
            if (n_Target_Nodes[typeIt][nbrIt] !=
                lttc->n_Target_Nodes[typeIt][nbrIt]) {
                fprintf(stderr, "n_Target_Nodes[] is incorrect\n");
                return (-1);
            }
        }
    }

    /* Verify lttc->target_Nodes[type][nbrs][nodes] are correct */
    for (typeIt = 0; typeIt < N_CELL_TYPES; typeIt++) {
        for (nbrIt = 0; nbrIt < lttc->n_Nbrs_Max; nbrIt++) {
            /* Make sure all elements after n_Target_Nodes[typeIt][nbrIt] are NULL */
            for (nodeIt = n_Target_Nodes[typeIt][nbrIt]; nodeIt < lttc->n_Nodes; nodeIt++) {
                if (lttc->target_Nodes[typeIt][nbrIt][nodeIt] != NULL) {
                    fprintf(stderr, "Not NULL cells exist in target_Nodes[%d][%lu] \
after n_Cells.\n", typeIt, nbrIt);
                    return (-1);
                }
            }
        }
    }

    for (typeIt = 0; typeIt < N_CELL_TYPES; typeIt++)
        free((void *) n_Target_Nodes[typeIt]);

    for (typeIt = NORMAL; typeIt < N_CELL_TYPES; typeIt++) {
        unsigned long int n_NbrsIt;
        for (n_NbrsIt = 0; n_NbrsIt <= lttc->n_Nbrs_Max; n_NbrsIt++) {
            for (nodeIt = 0; nodeIt < lttc->n_Target_Nodes[typeIt][n_NbrsIt]; nodeIt++) {
                if (lttc->target_Nodes[typeIt][n_NbrsIt][nodeIt]->type != TargetOf(typeIt))
                    return (-1);
            }
        }
    }


    return 0;
}/* End of VerifyState() */
