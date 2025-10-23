/**
 * @file packet_former.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for forming SCHP5 packets
 * @version 0.1
 * @date 2025-10-02
 */
# ifndef PACKET_FORMER_H
# define PACKET_FORMER_H

#include "packet.h"
#include "crcinterface.h"
#include <cstdint>

/**
 * @brief Forms packets according to the SCHP5 protocol for packet radios
 */
class PacketFormer {

    private:

    /**
    * @brief The address of the sender
    */
    const uint32_t senderAddress;

    /**
    * @brief The address of the reciever
    */
    const uint32_t recieverAddress;

    /**
     * @brief A pointer to the crc generator
     */
    crcutil_interface::CRC * crcGenny;

    /**
     * @brief A packet for storing information as output packets are formed
     */
    Packet tempPacket;


    public:

    /**
     * @brief Construct a new packet former
     */
    PacketFormer(
        uint8_t sender_address,
        uint8_t reciever_address);


    // Transmitter Packets

    /**
     * @brief Form a data packet
     * @param data A pointer to the input data buffer
     * @param sequence_number The sequence number of the current packet
     * @return A well formed SCHP5 data packet
     */
    Packet formDataPacket(std::vector<uint8_t> * data, uint8_t sequence_number, uint8_t data_length);

    /**
     * @brief Form a retransmitted packet
     * @param data A pointer to the input data buffer
     * @param sequence_number The sequence number of the current packet
     * @return A well formed SCHP5 retransmitted packet
     * 
     */
    Packet formRetransmitPacket(std::vector<uint8_t> * data, uint8_t sequence_num, uint8_t data_length);

    /**
     * @brief Form a retransmitted packet
     * @param data_packet A well formed SCHP5 data packet to retransmit
     * @return A well formed SCHP5 retransmitted packet
     * 
     */
    Packet formRetransmitPacket(Packet data_packet);

    /**
     * @brief Form a busy start packet
     * @return A well formed SCHP5 busy start packet
     */
    Packet formBusyStartPacket();

    /**
     * @brief Form a busy end packet
     * @return A well formed SCHP5 busy end packet
     */
    Packet formBusyEndPacket();

    /**
     * @brief Form a center frequency change packet
     * @param center_freq The center frequency to change to, as an IEEE 754 32-bit floating-point
     * @return A well formed SCHP5 center frequency change packet
     */
    Packet formCenterFrequencyPacket(const float center_freq);

    /**
     * @brief Form a modulation change packet
     * @param modulation_type The type of modulation to switch to, as specified in the SCHP5 protocol, usally log_2(QAM order - 1)
     * @return A well formed SCHP5 modulation change packet
     */
    Packet formModulationChangePacket(const uint8_t modulation_type);

    /**
     * @brief Form a transaction start packet
     * @param center_freq The center frequency to use, as an IEEE 754 32-bit floating-point
     * @param modulation_type The type of modulation to use, as specified in the SCHP5 protocol, usally log_2(QAM order - 1)
     * @return A well formed SCHP5 transaction start packet
     */
    Packet formTransactionStartPacket(const float center_freq, const uint8_t modulation_type);

    /**
     * @brief Form a transaction restart packet
     * @param center_freq The center frequency to use, as an IEEE 754 32-bit floating-point
     * @param modulation_type The type of modulation to use, as specified in the SCHP5 protocol, usally log_2(QAM order - 1)
     * @return A well formed SCHP5 transaction restart packet
     */
    Packet formTransactionRestartPacket(const float center_freq, const uint8_t modulation_type);

    /**
     * @brief Form a transaction transfer packet
     * @param center_freq The center frequency to use, as an IEEE 754 32-bit floating-point
     * @param modulation_type The type of modulation to use, as specified in the SCHP5 protocol, usally log_2(QAM order - 1)
     * @return A well formed SCHP5 transaction transfer packet
     */
    Packet formTransactionTransferPacket(const float center_freq, const uint8_t modulation_type);

    /**
     * @brief Form a transmitter transaction dropped packet
     * @return A well formed SCHP5 TX transmitter transaction dropped packet
     */
    Packet formTxTransactionDroppedPacket();

    /**
     * @brief Form a transaction end packet
     * @return A well formed SCHP5 transaction end packet
     */
    Packet formTransactionEndPacket();

    // Reciever Packets

    /**
     * @brief Form an acknowledge packet
     * @param error_count The current count of transmission errors
     * @param sequence_number The sequence number of the packet being acknowledged
     * @return A well formed SCHP5 acknowledge packet
     */
    Packet formAcknowledgePacket(uint16_t error_count, uint8_t sequence_number);

    /**
     * @brief Form a repeat request packet
     * @param sequence_number The sequence number of the packet being requested
     * @return A well formed SCHP5 repeat packet
     */
    Packet formRepeatPacket(uint8_t sequence_number);

    /**
     * @brief Form a busy accept packet
     * @return A well formed SCHP5 busy accept packet
     */
    Packet formBusyAcceptPacket();
    
    /**
     * @brief Form a busy ended packet
     * @return A well formed SCHP5 busy ended packet
     */
    Packet formBusyEndedPacket();

    /**
     * @brief Form a channel property change acknowledge packet
     * @return A well formed SCHP5 general channel acknowledge packet
     */
    Packet formGeneralChannelAcknowledge();

    /**
     * @brief Form a transaction start acknowledge packet
     * @return A well formed SCHP5 transaction start acknowledge packet
     */
    Packet formTransactionStartAcknowledge();

    /**
     * @brief Form a transfer accept packet
     * @return A well formed SCHP5 transfer accept packet
     */
    Packet formTransferAccept();

    /**
     * @brief Form a tranfer decline packet
     * @return A well formed SCHP5 transfer decline packet
     */
    Packet formTransferDecline();

    /**
     * @brief Form a reciever transaction dropped packet
     * @return A well formed SCHP5 RX transaction dropped packet
     */
    Packet formRxTransactionDroppedPacket();

    /**
     * @brief Form a transaction end acknowledge packet
     * @return A well formed SCHP5 transaction end acknowledge packet
     */
    Packet formTransactionEndAcknowledge();

};

#endif